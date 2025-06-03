rm(list = ls())

# Install and load mgcv if needed
if (!require(mgcv)) {
  install.packages("mgcv")
  library(mgcv)
}

# --------------------------
# Data Preparation
# --------------------------
# Load your new data
data <- read.csv("../data/8apr_final_data_wihtNewClimateVars.csv")
cat("Data dimensions:", dim(data), "\n")
head(data, 3)

# Create target variable: average CPB abundance
data$cpb_abundance <- (data$cpba_count + data$cpbl_count) / 2

# Convert croptype to binary (potato = 1 if crop code 43; else 0)
data$croptype <- ifelse(data$croptype == 43, 1, 0)

# Scale numeric predictors (adjust these based on your needs)
data$gdd                     <- scale(data$gdd)
data$cum_gdd                 <- scale(data$cum_gdd)
data$intensity               <- scale(data$intensity)
data$summer_avg_temp         <- scale(data$summer_avg_temp)
data$summer_avg_percip       <- scale(data$summer_avg_percip)
data$summer_heavy_rainfall_days <- scale(data$summer_heavy_rainfall_days)
data$summer_temp_variability <- scale(data$summer_temp_variability)
data$winter_avg_temp         <- scale(data$winter_extreme_cold_days)
data$winter_hottest_temp     <- scale(data$winter_hottest_temp)
data$winter_coldest_temp     <- scale(data$winter_coldest_temp)
data$winter_temp_variability <- scale(data$winter_temp_variability)
data$winter_warm_day_count   <- scale(data$winter_warm_day_count)
data$spring_frost_free_days  <- scale(data$spring_frost_free_days)
data$new_proportion          <- scale(data$new_proportion)

# Convert 'farm' and 'year' to factors for use as random effects
data$farm <- as.factor(data$farm)
data$year <- as.factor(data$year)

# --------------------------
# Fit GAM with Spatial Smoothing
# --------------------------
gam_model <- gam(cpb_abundance ~ s(lat, lng) + s(year, bs = "re") + s(farm, bs = "re") +
                   gdd + cum_gdd + croptype + intensity +
                   summer_avg_temp + summer_avg_percip + summer_heavy_rainfall_days + summer_temp_variability +
                   winter_extreme_cold_days + winter_hottest_temp + winter_coldest_temp + winter_temp_variability +
                   winter_warm_day_count + spring_frost_free_days + new_proportion,
                 data = data, family = gaussian())

# Print the model summary and AIC
summary(gam_model)
cat("AIC:", AIC(gam_model), "\n")

res <- residuals(gam_model)
mse <- mean(res^2)
print(mse)

# 1a. Fit the no‐spatial model
gam_nospatial <- gam(
  cpb_abundance ~ 
    s(year, bs = "re") + 
    s(farm, bs = "re") +
    gdd + cum_gdd + croptype + intensity +
    summer_avg_temp + summer_avg_percip + summer_heavy_rainfall_days + summer_temp_variability +
    winter_extreme_cold_days + winter_hottest_temp + winter_coldest_temp + winter_temp_variability +
    winter_warm_day_count + spring_frost_free_days + new_proportion,
  data = data, family = gaussian()
)

# 1b. Compare by AIC
AIC(gam_model, gam_nospatial)
# lower AIC = better

# —or— use a likelihood‐ratio test for nested GAMs
anova(gam_nospatial, gam_model, test = "Chisq")

# Plot the smooth effects (including the spatial effect)
plot(gam_model, scheme = 2, pages = 1, rug = TRUE, shade = TRUE)

