library(tidyverse)
library(lme4)
# Optional: for ICC calculation
if (!require(performance)) install.packages("performance")

library(performance)


# Input/output folder paths
RAW_DATA_DIR <- "../data/raw/"
INTERIM_DATA_DIR <- "../data/interim/"
PROCESSED_DATA_DIR <- "../data/processed/"
MODEL_DIR          <- "../models"

data <- read.csv(file.path(PROCESSED_DATA_DIR, 'final_data_for_modeling.csv'))
print(dim(data))

# Create target: average CPB abundance
data$cpb_abundance <- (data$cpba_count + data$cpbl_count) / 2

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

# find every column that is still a matrix, turn it into numeric
matrix_vars <- names(data)[sapply(data, is.matrix)]
for (v in matrix_vars) data[[v]] <- as.numeric(data[[v]])


### --------------------------------------------------------------------------------------------------------------------------------------
###  --------------------------------------------------  spaMM MODEL ---------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------------------------------

# --- Setup
if (!require(spaMM)) {
  install.packages("spaMM", repos = "https://pbil.univ-lyon1.fr/CRAN/")
  library(spaMM)
}

spamm_SRF <- fitme(
  cpb_abundance ~ cum_gdd + wei_intensity + summer_heavy_rainfall_days + winter_coldest_temp +
    Matern(1 | lat + lng) + (1 | year) + (1 | farm),
  data = data, family = gaussian()
)

summary(spamm_SRF)

AIC_SRF          <- AIC(spamm_SRF)
res_SRF          <- residuals(spamm_SRF)
MSE_SRF          <- mean(res_SRF^2)
var_fixed_SRF    <- var(predict(spamm_SRF))
var_spatial_SRF  <- unname(spamm_SRF$lambda["lat + lng"])
var_year_SRF     <- unname(spamm_SRF$lambda["year"])
var_farm_SRF     <- unname(spamm_SRF$lambda["farm"])
var_resid_SRF    <- unname(spamm_SRF$phi)
R2m_SRF          <- var_fixed_SRF / (var_fixed_SRF + var_spatial_SRF + var_year_SRF + var_farm_SRF + var_resid_SRF)
R2c_SRF          <- (var_fixed_SRF + var_spatial_SRF + var_year_SRF + var_farm_SRF) /
  (var_fixed_SRF + var_spatial_SRF + var_year_SRF + var_farm_SRF + var_resid_SRF)

cat("\n--- Model spamm_SRF (Fixed + Random + Spatial) ---\n")

cat("AIC:           ", round(AIC_SRF, 2), "\n")
cat("MSE:           ", round(MSE_SRF, 3), "\n")
cat("R² (marginal): ", round(R2m_SRF, 3), "\n")
cat("R² (conditional): ", round(R2c_SRF, 3), "\n")


spamm_red <- fitme(
  cpb_abundance ~ cum_gdd + wei_intensity + summer_heavy_rainfall_days + winter_coldest_temp +
    Matern(1 | lat + lng) + (1 | year),
  data = data, family = gaussian()
)
summary(spamm_red)

saveRDS(spamm_red, file.path(MODEL_DIR, "spamm_red_noFarm.rds"))

print(AIC(spamm_red))
# --- MSE ---
res_red <- residuals(spamm_red)
MSE_red <- mean(res_red^2)

# --- Variance components ---
var_fixed_red   <- var(predict(spamm_red))
var_spatial_red <- unname(spamm_red$lambda["lat + lng"])
var_year_red    <- unname(spamm_red$lambda["year"])
var_resid_red   <- unname(spamm_red$phi)

# --- R² ---
R2m_red <- var_fixed_red / (var_fixed_red + var_spatial_red + var_year_red + var_resid_red)
R2c_red <- (var_fixed_red + var_spatial_red + var_year_red) /
  (var_fixed_red + var_spatial_red + var_year_red + var_resid_red)

# --- Print ---
cat("\n--- Model spamm_red (No farm RE) ---\n")
cat("AIC:           ", round(AIC(spamm_red), 2), "\n")
cat("MSE:           ", round(MSE_red, 3), "\n")
cat("R² (marginal): ", round(R2m_red, 3), "\n")
cat("R² (conditional): ", round(R2c_red, 3), "\n")
