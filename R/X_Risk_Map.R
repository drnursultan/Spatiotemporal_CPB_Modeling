# ===== Packages =====
libs <- c("dplyr","readr","sf","ggplot2","viridis","spaMM","tigris")
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

# ===== 1) Load saved model + training stats for scaling =====
spamm_red <- readRDS(file.path(MODEL_DIR, "spamm_red_noFarm.rds"))

train <- read_csv(file.path(PROCESSED_DATA_DIR, "final_data_for_modeling.csv"),
                  show_col_types = FALSE) %>%
  mutate(
    cpb_abundance = (cpba_count + cpbl_count)/2,
    year = factor(year)
  )

# Predictors used in the saved model
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

# ===== 2) Load grid attributes (already assembled) =====
# Expected columns: grid10_id, lon, lat, cum_gdd, wei_intensity, summer_heavy_rainfall_days, winter_coldest_temp
grid_attr <- read_csv(file.path(PROCESSED_DATA_DIR, "grid10km_final_risk_attributes.csv"),
                      show_col_types = FALSE)

need_cols <- c("grid10_id","lon","lat", pred_vars)
miss <- setdiff(need_cols, names(grid_attr))
if (length(miss)) stop("Missing columns in grid file: ", paste(miss, collapse=", "))

grid_pred <- grid_attr %>%
  transmute(
    grid10_id,
    lng = lon,     # spaMM expects 'lng' + 'lat' for Matern()
    lat = lat,
    across(all_of(pred_vars))
  )

# Fill any NAs with training means, then scale with the training recipe
for (v in pred_vars) {
  if (anyNA(grid_pred[[v]])) grid_pred[[v]][is.na(grid_pred[[v]])] <- scalers[[v]]$mean
  grid_pred[[v]] <- scale_var(grid_pred[[v]], v)
}

# Year factor for prediction (match training levels)
grid_pred$year <- factor("2021", levels = levels(train$year))

# ===== 3) Predict =====
pred <- predict(spamm_red, newdata = grid_pred, type = "response")
grid_pred$pred_abund <- as.numeric(pred)

# Quick diagnostics
cat("Non-finite predictions: ", sum(!is.finite(grid_pred$pred_abund)), "\n")
print(summary(grid_pred$pred_abund))
print(quantile(grid_pred$pred_abund, probs = c(0, .01, .5, .99, 1), na.rm = TRUE))

# ===== 4) Build RiskScore safely =====
finite <- is.finite(grid_pred$pred_abund)

# Min–max scale using only finite values
minv <- min(grid_pred$pred_abund[finite], na.rm = TRUE)
maxv <- max(grid_pred$pred_abund[finite], na.rm = TRUE)
span <- maxv - minv
if (!is.finite(span) || span == 0) {
  warning("Predictions have zero or non-finite range; RiskScore set to NA.")
  grid_pred$RiskScore <- NA_real_
} else {
  grid_pred$RiskScore <- NA_real_
  grid_pred$RiskScore[finite] <- (grid_pred$pred_abund[finite] - minv) / span
}

# Optional: quantile stretch (better visual contrast, robust to outliers)
q <- quantile(grid_pred$pred_abund[finite], probs = c(0.02, 0.98), na.rm = TRUE)
grid_pred$RiskScore_q <- pmin(1, pmax(0, (grid_pred$pred_abund - q[1]) / (q[2] - q[1])))

# ===== 5) Join to 10-km polygons =====
g10_sf <- st_read(file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.gpkg"), quiet = TRUE)

g10_pred <- g10_sf %>%
  left_join(grid_pred %>% select(grid10_id, pred_abund, RiskScore, RiskScore_q), by = "grid10_id")

# ===== 6) Plot WI map =====
wi_sf <- states(cb = TRUE) %>%
  dplyr::filter(STUSPS == "WI") %>%
  st_transform(st_crs(g10_pred))

p_wi_10km <- ggplot() +
  geom_sf(data = wi_sf, fill = "grey95", color = "white", linewidth = 0.2) +
  geom_sf(data = g10_pred, aes(fill = RiskScore_q), color = NA) +  # use RiskScore_q for better contrast
  scale_fill_viridis(name = "Risk (0–1)", limits = c(0,1), na.value = "grey85") +
  labs(title = "CPB Risk (2021) – Wisconsin (10 km grid)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

print(p_wi_10km)
ggsave(file.path(PLOT_DIR, "CPB_Risk_WI_10km_2021_noFarm.png"),
       p_wi_10km, width = 9, height = 7, dpi = 300)

# ===== 7) Save predictions table (optional) =====
pred_out <- grid_pred %>%
  select(grid10_id, lng, lat, pred_abund, RiskScore, RiskScore_q)

write_csv(pred_out, file.path(PROCESSED_DATA_DIR, "grid10km_predictions_2021_noFarm.csv"))

# Peek
pred_out %>% head() %>% print()
