library(tidyverse)
library(lme4)
# Optional: for ICC calculation
if (!require(performance)) install.packages("performance")

library(performance)


# Input/output folder paths
RAW_DATA_DIR <- "../data/raw/"
INTERIM_DATA_DIR <- "../data/interim/"
PROCESSED_DATA_DIR <- "../data/processed/"

data <- read.csv(file.path(PROCESSED_DATA_DIR, 'final_data_for_modeling.csv'))
print(dim(data))

# Create target: average CPB abundance
data$cpb_abundance <- (data$cpba_count + data$cpbl_count) / 2
data$croptype <- ifelse(data$croptype == 43, 1, 0)


# --- Scale numeric predictors
numeric_vars <- c(
  "gdd", "cum_gdd", "wei_intensity", "wei_prop",
  "summer_avg_temp", "summer_avg_percip", "summer_heavy_rainfall_days",
  "summer_hottest_temp","summer_temp_variability", "winter_coldest_temp",
  "winter_extreme_cold_days",
  "winter_warm_day_count"
)

# Scale each numeric variable if it exists in data
for (v in numeric_vars) {
  if (v %in% names(data)) data[[v]] <- scale(data[[v]])
}

nrow(data[data['croptype']==1,])

# find every column that is still a matrix, turn it into numeric
matrix_vars <- names(data)[sapply(data, is.matrix)]
for (v in matrix_vars) data[[v]] <- as.numeric(data[[v]])

### --------------------------------------------------------------------------------------------------------------------------------------
###  --------------------------------------------------  LMER MODEL ----------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------------------------------

# (a) Null vs. farm random effect
model_null  <- lm(cpb_abundance ~ 1, data = data)
model_farm  <- lmer(cpb_abundance ~ 1 + (1 | farm), data = data, REML = TRUE)
anova(model_farm,model_null)

# (b) Farm vs. farm+year random effects
model_farm_year <- lmer(cpb_abundance ~ 1 + (1 | farm) + (1 | year), data = data, REML = TRUE)
anova(model_farm, model_farm_year)

# --- Full mixed-effects model
full_model <- lmer(
  cpb_abundance ~ gdd + cum_gdd + croptype + wei_intensity + wei_prop +
    summer_avg_temp + summer_avg_percip + summer_heavy_rainfall_days +
    summer_temp_variability + summer_hottest_temp + winter_coldest_temp + winter_extreme_cold_days +
    winter_warm_day_count +
    (1 | farm) + (1 | year),
  data = data, REML = TRUE
)

summary(full_model)

### ------------------------------------------------------  PLOTS  ----------------------------------------------------------------------

#install.packages("insight")       # ensures you have ≥ 1.2.0
library(sjPlot)                   # should now load insight ≥ 1.2.0
library(performance)              # for check_model()

# basic caterpillar plot of all random effects
plot_model(full_model, type = "re", 
           sort.est = TRUE,        # orders by size
           show.values = TRUE,     # prints the estimates
           value.offset = .3)      # nudges the labels


# 1. Grab your diagnostics
plot(check_model(full_model, panel = FALSE))

# extract the farm random intercepts
farm_intercepts <- ranef(full_model)$farm[,1]

# simple horizontal boxplot
boxplot(
  farm_intercepts,
  horizontal = TRUE,
  main       = "Distribution of Farm Random Intercepts",
  xlab       = "Random intercept"
)


### --------------------------------------------------------------------------------------------------------------------------------------
###  --------------------------------------------------  GAMM MODEL ----------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------------------------------

# --- Setup
if (!require(mgcv)) {
  install.packages("mgcv")
  library(mgcv)
}

data$farm <- as.factor(data$farm)
data$year <- as.factor(data$year)

gam_model <- gam(cpb_abundance ~ s(lat, lng) + s(year, bs = "re") + s(farm, bs = "re") +
                   gdd + cum_gdd + croptype + wei_intensity + wei_prop +
                   summer_avg_temp + summer_avg_percip + summer_heavy_rainfall_days +
                   summer_temp_variability + summer_hottest_temp + winter_coldest_temp + winter_extreme_cold_days +
                   winter_warm_day_count,
                 data = data, family = gaussian())

# Print the model summary and AIC
summary(gam_model)

### ------------------------------------------------------  PLOTS  ----------------------------------------------------------------------


# ──────────────────────────────────────────────────────────────────────────
# 0. Load packages (install if needed)
# ──────────────────────────────────────────────────────────────────────────
if (!require(mgcv))      install.packages("mgcv")  
library(mgcv)

if (!require(ggplot2))   install.packages("ggplot2")
library(ggplot2)

# ──────────────────────────────────────────────────────────────────────────
# 1. Smooth‐term plots
#───────────────────────────────────────────────────────────────────────────
# This will plot each smooth (including your spatial s(lat,lng) and the two
# random‐effect smooths) in separate panels.  `shade=TRUE` draws ±2·SE bands.
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))  
plot(gam_model, select = 1, shade = TRUE, main = "s(lat, lng)")  
plot(gam_model, select = 2, shade = TRUE, main = "s(year) [RE]")  
plot(gam_model, select = 3, shade = TRUE, main = "s(farm) [RE]")  
par(mfrow = c(1, 1))  # reset

# ──────────────────────────────────────────────────────────────────────────
# 2. Fitted vs. Observed
#───────────────────────────────────────────────────────────────────────────
# Extract fitted values and plot against the true response with a 1:1 line.
fitted_vals <- predict(gam_model, type = "response")
obs        <- gam_model$y

ggplot(data.frame(obs, fitted_vals), aes(x = obs, y = fitted_vals)) +
  geom_point(alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    x     = "Observed cpb_abundance",
    y     = "Fitted values",
    title = "Fitted vs. Observed"
  ) +
  theme_minimal()

# ──────────────────────────────────────────────────────────────────────────
# 3. Residual diagnostics
#───────────────────────────────────────────────────────────────────────────
resid_vals <- residuals(gam_model, type = "response")
fitted_vals

# arrange 2×2 grid
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# 3a. Residuals vs Fitted
plot(fitted_vals, resid_vals,
     xlab = "Fitted", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

# 3b. QQ‐plot of residuals
qqnorm(resid_vals, main = "QQ‐plot")
qqline(resid_vals)

# 3c. Histogram of residuals
hist(resid_vals,
     main = "Histogram of Residuals", xlab = "Residual")

# 3d. Scale‐location (sqrt of |resid|)
plot(fitted_vals, sqrt(abs(resid_vals)),
     xlab = "Fitted", ylab = "√|Residual|",
     main = "Scale–Location")
par(mfrow = c(1, 1))  # reset

# ──────────────────────────────────────────────────────────────────────────
# 4. Heatmap of spatial smooth s(lat, lng)
#───────────────────────────────────────────────────────────────────────────
# Use mgcv::vis.gam to draw a contour (heatmap) of the 2D smooth.
# `type="terms"` plots just the spatial term, holding others at their
# reference/mean values.
# 4. Heatmap of spatial smooth s(lat, lng) on the response scale
vis.gam(
  gam_model,
  view      = c("lat", "lng"),    # variables on the two axes
  plot.type = "contour",          # filled contour; use "persp" for 3D
  type      = "response",         # "link" or "response"
  se        = TRUE,               # add ±SE contours
  too.far   = 0.1,                # drop points outside convex hull
  select    = 1,                  # 1st smooth term is s(lat,lng)
  main      = "Spatial smooth: s(lat, lng)",
  xlab      = "Latitude",
  ylab      = "Longitude"
)


### --------------------------------------------------------------------------------------------------------------------------------------
###  --------------------------------------------------  spaMM MODEL ---------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------------------------------

# --- Setup
if (!require(spaMM)) {
  install.packages("spaMM", repos = "https://pbil.univ-lyon1.fr/CRAN/")
  library(spaMM)
}

spamm_SRF <- fitme(
  cpb_abundance ~
    gdd + cum_gdd + croptype + wei_intensity + wei_prop +
    summer_avg_temp + summer_avg_percip + summer_heavy_rainfall_days +
    summer_temp_variability + summer_hottest_temp + winter_coldest_temp + winter_extreme_cold_days +
    winter_warm_day_count +
    Matern(1 | lat + lng) + (1 | year) + (1 | farm),
  data = data, family = gaussian()
)

summary(spamm_SRF)

### ------------------------------------------------------  PLOTS  ----------------------------------------------------------------------

# ─────────────────────────────────────────────────────────────────────────────
# 0. Packages (install if needed)
# ─────────────────────────────────────────────────────────────────────────────
pkgs <- c("spaMM","ggplot2","akima","viridis","gstat","sp","forcats")
for (p in pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    install.packages(p, repos = "https://cran.rstudio.com")
  }
  library(p, character.only=TRUE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. Spatial Random-Effect Map (kriging-style)
# ─────────────────────────────────────────────────────────────────────────────

# 1. pull out the raw list of ranefs
rL <- ranef(spamm_SRF)

# 2. extract the three raw vectors
mat_re   <- rL[["Matern(1 | lat + lng)"]]   # length = # unique coords
year_re  <- rL[["( 1 | year )"]]            # length = # levels(year)
farm_re  <- rL[["( 1 | farm )"]]            # length = # levels(farm)

# 3. rebuild the unique‐coords data.frame in the SAME order spaMM used:
uniq_coords <- unique(data[c("lat","lng")])  # first‐occurrence order
stopifnot(nrow(uniq_coords)==length(mat_re))

# 4. attach the Matern BLUPs to those coords
uniq_coords$sp_re <- mat_re

# 5. now look up each observation’s spatial effect via a join‐by:
key_obs  <- paste(data$lat, data$lng, sep=";")
key_uniq <- paste(uniq_coords$lat, uniq_coords$lng, sep=";")
spat_re  <- uniq_coords$sp_re[ match(key_obs, key_uniq) ]
stopifnot(length(spat_re)==nrow(data))  # sanity check

# 6. build the coords_df for plotting
coords_df <- data.frame(
  lng   = data$lng,
  lat   = data$lat,
  sp_re = spat_re
)

# 7. krige‐style heatmap
library(akima)
interp_grid <- interp(
  x         = coords_df$lng,
  y         = coords_df$lat,
  z         = coords_df$sp_re,
  duplicate = "mean", nx = 200, ny = 200
)
grid_df <- expand.grid(lng = interp_grid$x, lat = interp_grid$y)
grid_df$sp_re <- as.vector(interp_grid$z)

library(ggplot2); library(viridis)
ggplot(grid_df, aes(lng, lat, fill=sp_re)) +
  geom_raster() +
  scale_fill_viridis_c(option="D") +
  coord_quickmap() +
  labs(
    title = "Spatial Random-Effect Map (Matern RF)",
    x = "Longitude", y = "Latitude", fill = "Spatial RE"
  ) +
  theme_minimal()

head(mat_re)         # first 6 spatial BLUPs
head(uniq_coords)    # the corresponding lat/lng


# ---------------------------------------------------------------------------------------------------

library(ggplot2)
library(forcats)

# 1. Farm intercepts
farm_re <- ranef(spamm_SRF)[["( 1 | farm )"]]      # a named numeric vector
farm_df <- data.frame(
  farm      = factor(names(farm_re), levels = names(farm_re)),
  intercept = as.numeric(farm_re)
)
farm_df$farm <- fct_reorder(farm_df$farm, farm_df$intercept)

ggplot(farm_df, aes(x = farm, y = intercept)) +
  geom_point() +
  coord_flip() +
  labs(
    title = "Farm Random Effects (Intercepts)",
    x     = "Farm",
    y     = "Deviation from overall mean"
  ) +
  theme_minimal()


# 2. Year intercepts
year_re <- ranef(spamm_SRF)[["( 1 | year )"]]      # another named vector
year_df <- data.frame(
  year      = factor(names(year_re), levels = names(year_re)),
  intercept = as.numeric(year_re)
)
year_df$year <- fct_reorder(year_df$year, year_df$intercept)

ggplot(year_df, aes(x = year, y = intercept)) +
  geom_point() +
  coord_flip() +
  labs(
    title = "Year Random Effects (Intercepts)",
    x     = "Year",
    y     = "Deviation from overall mean"
  ) +
  theme_minimal()

# ---------------------------------------------------------------------------------------------------


library(sp); library(gstat)

resid_df <- data.frame(
  lng   = data$lng,
  lat   = data$lat,
  resid = residuals(spamm_SRF, type = "response")
)
coordinates(resid_df) <- ~ lng + lat

vg_emp <- variogram(resid ~ 1, resid_df)
plot(vg_emp$dist, vg_emp$gamma,
     pch   = 16,
     xlab  = "Distance",
     ylab  = "Semivariance",
     main  = "Empirical Variogram of Residuals")
