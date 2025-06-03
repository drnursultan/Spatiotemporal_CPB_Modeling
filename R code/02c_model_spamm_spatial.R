rm(list = ls())

# Install and load required packages (if not already installed)
if (!require(spaMM)) {
  install.packages("spaMM", repos = "https://pbil.univ-lyon1.fr/CRAN/")
  library(spaMM)
}
if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

# --------------------------
# 1. Data Preparation
# --------------------------
data <- read.csv("../output_new_data/20may_final_data.csv")

# Create target: average CPB abundance
data$cpb_abundance <- (data$cpba_count + data$cpbl_count) / 2

# Scale numeric predictors
data$gdd                        <- scale(data$gdd)
data$cum_gdd                    <- scale(data$cum_gdd)
data$intensity                  <- scale(data$intensity)
data$summer_avg_temp            <- scale(data$summer_avg_temp)
data$summer_avg_percip          <- scale(data$summer_avg_percip)
data$summer_heavy_rainfall_days <- scale(data$summer_heavy_rainfall_days)
data$summer_temp_variability    <- scale(data$summer_temp_variability)
data$winter_avg_temp            <- scale(data$winter_extreme_cold_days)
data$winter_hottest_temp        <- scale(data$winter_hottest_temp)
data$winter_coldest_temp        <- scale(data$winter_coldest_temp)
data$winter_temp_variability    <- scale(data$winter_temp_variability)
data$winter_warm_day_count      <- scale(data$winter_warm_day_count)
data$spring_frost_free_days     <- scale(data$spring_frost_free_days)
data$new_proportion             <- scale(data$new_proportion)

# Convert 'farm' and 'year' to factors for random effects
data$farm <- as.factor(data$farm)
data$year <- as.factor(data$year)

# --------------------------
# 2. Model spamm_RF: Fixed + Random (year + farm)
# --------------------------

spamm_RF <- fitme(
  cpb_abundance ~
    gdd + cum_gdd + croptype + intensity +
    summer_avg_temp + summer_avg_percip +
    summer_heavy_rainfall_days + summer_temp_variability +
    winter_extreme_cold_days + winter_hottest_temp + winter_coldest_temp + winter_temp_variability +
    winter_warm_day_count + spring_frost_free_days + new_proportion +
    (1 | year) + (1 | farm),
  data = data, family = gaussian()
)

# 2a. Metrics for spamm_RF
AIC_RF         <- AIC(spamm_RF)
res_RF         <- residuals(spamm_RF)
MSE_RF         <- mean(res_RF^2)
var_fixed_RF   <- var(predict(spamm_RF))
var_year_RF    <- unname(spamm_RF$lambda["year"])
var_farm_RF    <- unname(spamm_RF$lambda["farm"])
var_resid_RF   <- unname(spamm_RF$phi)
R2m_RF         <- var_fixed_RF / (var_fixed_RF + var_year_RF + var_farm_RF + var_resid_RF)
R2c_RF         <- (var_fixed_RF + var_year_RF + var_farm_RF) /
  (var_fixed_RF + var_year_RF + var_farm_RF + var_resid_RF)

cat("\n--- Model spamm_RF (Fixed + Random) ---\n")
cat("AIC:           ", round(AIC_RF, 2), "\n")
cat("MSE:           ", round(MSE_RF, 3), "\n")
cat("R² (marginal): ", round(R2m_RF, 3), "\n")
cat("R² (conditional): ", round(R2c_RF, 3), "\n")


# --------------------------
# 3. Model spamm_SRF: Fixed + Random + Spatial
# --------------------------

spamm_SRF <- fitme(
  cpb_abundance ~
    gdd + cum_gdd + croptype + intensity +
    summer_avg_temp + summer_avg_percip +
    summer_heavy_rainfall_days + summer_temp_variability +
    winter_extreme_cold_days + winter_hottest_temp + winter_coldest_temp + winter_temp_variability +
    winter_warm_day_count + spring_frost_free_days + new_proportion +
    Matern(1 | lat + lng) + (1 | year) + (1 | farm),
  data = data, family = gaussian()
)

# 3a. Metrics for spamm_SRF
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


# 1. Summaries
summary(spamm_RF)
summary(spamm_SRF)


# 1. Extract the “conditional AIC” (2nd element) from each AIC vector
condAIC_RF  <- AIC_RF[2]
condAIC_SRF <- AIC_SRF[2]

# 2. Print them
cat(
  "Conditional AIC (spamm_RF):  ", round(condAIC_RF, 2), "\n",
  "Conditional AIC (spamm_SRF): ", round(condAIC_SRF, 2), "\n"
)

# 3. Run the bootstrap LRT and capture its p‑value
lrt     <- anova(spamm_RF, spamm_SRF, boot.repl = 200)
p_boot  <- lrt$`Pr(>Chisq)`[2]

cat("Bootstrap LRT p‑value for adding spatial term:", signif(p_boot, 3), "\n")
