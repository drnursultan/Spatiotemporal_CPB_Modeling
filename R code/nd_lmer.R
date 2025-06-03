rm(list = ls())

install.packages("lme4")
install.packages("performance")
install.packages("ggplot2")

library(lme4)
library(performance)
library(ggplot2)

# --------------------------
# 1. Load New Data and Prepare Variables
# --------------------------
data <- read.csv("../data/8apr_final_data_wihtNewClimateVars.csv")
print(dim(data))
print(head(data, 3))

# Create target: average CPB abundance
data$cpb_abundance <- (data$cpba_count + data$cpbl_count) / 2

# Convert croptype to binary (1 if potato, crop code 43; 0 otherwise)
data$croptype <- ifelse(data$croptype == 43, 1, 0)

# Scale numeric predictors (you may adjust which variables you scale)
data$gdd                     <- scale(data$gdd)
data$cum_gdd                 <- scale(data$cum_gdd)
data$intensity               <- scale(data$intensity)
#data$potato_proportion       <- scale(data$potato_proportion)
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

# --------------------------
# 2. Random Effects Testing
# --------------------------
# Since the new dataset does not have 'grower', we compare models with random effects of farm and year.

# (a) Test the importance of the random effect of farm
model_farm <- lmer(cpb_abundance ~ 1 + (1 | farm), data = data, REML = TRUE)
model_null <- lm(cpb_abundance ~ 1, data = data)
anova(model_farm, model_null)

# (b) Test if adding year as a random effect improves the model
model_farm_year <- lmer(cpb_abundance ~ 1 + (1 | farm) + (1 | year), data = data, REML = TRUE)
anova(model_farm, model_farm_year)

# --------------------------
# 3. Build the Full Mixed-Effects Model (Fixed Effects + Random Effects)
# --------------------------
full_model <- lmer(cpb_abundance ~ gdd + cum_gdd + croptype + intensity +
                     summer_avg_temp + summer_avg_percip + summer_heavy_rainfall_days + summer_temp_variability +
                     winter_extreme_cold_days + winter_hottest_temp + winter_coldest_temp + winter_temp_variability +
                     winter_warm_day_count + spring_frost_free_days + new_proportion +
                     (1 | farm) + (1 | year), 
                   data = data, REML = TRUE)

summary(full_model)

# Compute ICC to see how much variance is explained by the random effects (farm and year)
icc_val <- icc(full_model)
print(icc_val)

# Display model AIC and Mean Squared Error of residuals
print(AIC(full_model))
mse <- mean(residuals(full_model)^2)
print(mse)

# --------------------------
# 4. Diagnostic Plots
# --------------------------
# Residual plot
plot(residuals(full_model), main = "Residual Plot", ylab = "Residuals", xlab = "Index")
# QQ plot for residuals
qqnorm(residuals(full_model), main = "QQ Plot of Residuals")
qqline(residuals(full_model))


# testing for spatial corr:

# 1. Compute residuals from your fitted mixed model
resid_vals <- residuals(full_model)

# 2. Build a spatial neighbors object from your field coordinates
install.packages("spdep")
library(spdep)
coords <- cbind(data$lng, data$lat)

# e.g. k = 4 nearest neighbors
knn_nb  <- knearneigh(coords, k = 4)
nb_list <- knn2nb(knn_nb)

# row‐standardized weights
lw <- nb2listw(nb_list, style = "W")

# 3. Moran’s I test for global spatial autocorrelation
moran_res <- moran.test(resid_vals, lw)
print(moran_res)
