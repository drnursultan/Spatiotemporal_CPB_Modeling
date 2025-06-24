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
