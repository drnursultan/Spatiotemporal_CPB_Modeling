# ===== Packages =====
libs <- c("dplyr","readr","sf","ggplot2","viridis","spaMM","tigris","purrr","stringr")
invisible(lapply(libs, function(p){
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))
options(tigris_use_cache = TRUE, tigris_class = "sf")

# ===== Paths =====
INTERIM_DATA_DIR   <- "../data/interim/"
PROCESSED_DATA_DIR <- "../data/processed/"
MODEL_DIR          <- "../models"
PLOT_DIR           <- "../plots/maps"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ===== 1) Load model + training (for scalers) =====
spamm_red <- readRDS(file.path(MODEL_DIR, "spamm_red_noFarm.rds"))

train <- read_csv(file.path(PROCESSED_DATA_DIR, "final_data_for_modeling.csv"),
                  show_col_types = FALSE) %>%
  mutate(
    cpb_abundance = (cpba_count + cpbl_count)/2,
    year = factor(year)
  )

# Predictors used in the model
pred_vars <- c("cum_gdd","wei_intensity","summer_heavy_rainfall_days","winter_coldest_temp")

# Training means/SDs for scaling
scalers <- lapply(pred_vars, function(v){
  stopifnot(v %in% names(train))
  m <- mean(train[[v]], na.rm = TRUE); s <- sd(train[[v]], na.rm = TRUE)
  list(var=v, mean=m, sd=s)
})
names(scalers) <- pred_vars

scale_var <- function(x, v){
  m <- scalers[[v]]$mean; s <- scalers[[v]]$sd
  if (!is.finite(s) || s == 0) x - m else (x - m)/s
}

# ===== 2) Base grid attributes =====
grid_attr <- read_csv(file.path(PROCESSED_DATA_DIR, "grid10km_final_risk_attributes.csv"),
                      show_col_types = FALSE)

need_cols <- c("grid10_id","lon","lat", pred_vars)
miss <- setdiff(need_cols, names(grid_attr))
if (length(miss)) stop("Missing columns in grid file: ", paste(miss, collapse=", "))

grid_base <- grid_attr %>%
  transmute(
    grid10_id,
    lng = lon,
    lat = lat,
    across(all_of(pred_vars))
  )

# ===== 3) Scenario definitions =====
# Each function takes a df and returns a modified df (unscaled).
scenarios <- list(
  baseline = function(df) df,
  high_intensity_1p5x = function(df) {
    df$wei_intensity <- pmin(1, df$wei_intensity * 1.5)
    df
  },
  low_intensity_0p7x = function(df) {
    df$wei_intensity <- pmax(0, df$wei_intensity * 0.7)
    df
  },
  warmer_winter_plus2C = function(df) {
    df$winter_coldest_temp <- df$winter_coldest_temp + 2
    df
  },
  colder_winter_minus2C = function(df) {
    df$winter_coldest_temp <- df$winter_coldest_temp - 2
    df
  },
  wetter_summer_heavyRain_plus2 = function(df) {
    df$summer_heavy_rainfall_days <- pmax(0, df$summer_heavy_rainfall_days + 2)
    df
  },
  drier_summer_heavyRain_minus2 = function(df) {
    df$summer_heavy_rainfall_days <- pmax(0, df$summer_heavy_rainfall_days - 2)
    df
  },
  warmer_growing_cumgdd_plus10pct = function(df) {
    df$cum_gdd <- df$cum_gdd * 1.10
    df
  },
  cooler_growing_cumgdd_minus10pct = function(df) {
    df$cum_gdd <- df$cum_gdd * 0.90
    df
  }
)

# ===== 4) Predict helper =====
predict_grid <- function(df_unscaled, year_level = "2021") {
  df <- df_unscaled
  
  # Fill NAs with training means, then scale
  for (v in pred_vars) {
    if (anyNA(df[[v]])) df[[v]][is.na(df[[v]])] <- scalers[[v]]$mean
    df[[v]] <- scale_var(df[[v]], v)
  }
  df$year <- factor(year_level, levels = levels(train$year))
  
  # Predict
  pred <- predict(spamm_red, newdata = df, type = "response")
  df$pred_abund <- as.numeric(pred)
  
  # Risk scores
  finite <- is.finite(df$pred_abund)
  minv <- min(df$pred_abund[finite], na.rm = TRUE)
  maxv <- max(df$pred_abund[finite], na.rm = TRUE)
  span <- maxv - minv
  df$RiskScore <- NA_real_
  if (is.finite(span) && span > 0) {
    df$RiskScore[finite] <- (df$pred_abund[finite] - minv) / span
  }
  
  q <- quantile(df$pred_abund[finite], probs = c(0.02, 0.98), na.rm = TRUE)
  df$RiskScore_q <- pmin(1, pmax(0, (df$pred_abund - q[1]) / (q[2] - q[1])))
  
  df
}

# ===== 5) Wisconsin layer for plots =====
g10_sf <- st_read(file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.gpkg"), quiet = TRUE)
wi_sf  <- states(cb = TRUE) %>%
  dplyr::filter(STUSPS == "WI") %>%
  st_transform(st_crs(g10_sf))

# ===== 6) Run scenarios, map, and save =====
results <- list()

for (nm in names(scenarios)) {
  cat("Running scenario:", nm, "\n")
  df_mod <- scenarios[[nm]](grid_base)
  pred_df <- predict_grid(df_mod, year_level = "2021") %>%
    select(grid10_id, lng, lat, all_of(pred_vars), pred_abund, RiskScore, RiskScore_q)
  
  # Join to polygons
  map_sf <- g10_sf %>%
    left_join(pred_df %>% select(grid10_id, pred_abund, RiskScore_q), by = "grid10_id")
  
  # Plot (using quantile-stretched score)
  p <- ggplot() +
    geom_sf(data = wi_sf, fill = "grey95", color = "white", linewidth = 0.2) +
    geom_sf(data = map_sf, aes(fill = RiskScore_q), color = NA) +
    scale_fill_viridis(name = "Risk (0–1)", limits = c(0,1), na.value = "grey85") +
    labs(title = paste0("CPB Risk (2021) – 10 km grid – ", nm)) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  
  ggsave(file.path(PLOT_DIR, paste0("CPB_Risk_WI_10km_2021_", nm, ".png")),
         p, width = 9, height = 7, dpi = 300)
  
  # Save CSV per scenario
  write_csv(pred_df, file.path(PROCESSED_DATA_DIR, paste0("grid10km_predictions_2021_", nm, ".csv")))
  
  results[[nm]] <- pred_df
}

# ===== 7) Quick comparison vs baseline =====
baseline <- results$baseline %>% select(grid10_id, RiskScore_q) %>%
  rename(RiskScore_q_baseline = RiskScore_q)

summary_tbl <- bind_rows(lapply(names(results), function(nm){
  df <- results[[nm]] %>% select(grid10_id, pred_abund, RiskScore_q) %>% mutate(scenario = nm)
})) %>%
  left_join(baseline, by = "grid10_id") %>%
  mutate(delta_q = RiskScore_q - RiskScore_q_baseline) %>%
  group_by(scenario) %>%
  summarise(
    n = n(),
    mean_pred = mean(pred_abund, na.rm = TRUE),
    median_pred = median(pred_abund, na.rm = TRUE),
    mean_risk_q = mean(RiskScore_q, na.rm = TRUE),
    median_risk_q = median(RiskScore_q, na.rm = TRUE),
    mean_delta_q_vs_baseline = mean(delta_q, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(scenario, names(scenarios)))  # keep scenario order

print(summary_tbl)
write_csv(summary_tbl, file.path(PROCESSED_DATA_DIR, "scenario_summary_2021.csv"))

