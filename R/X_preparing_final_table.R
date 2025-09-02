# make_final_risk_table.R
# Build the final attribute table for risk map (join weather + CDL intensity/prop)

library(dplyr)
library(readr)
library(sf)

# --- Paths ---
INTERIM_DATA_DIR   <- "../data/interim/"
PROCESSED_DATA_DIR <- "../data/processed/"

# Input files (change names if yours differ)
gpkg_grid_path   <- file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.gpkg")
weather_csv_path <- file.path(INTERIM_DATA_DIR, "grid10km_climate_2021.csv")      # <- your saved weather table (with gdd/cum_gdd/etc.)
cdl_csv_path     <- file.path(INTERIM_DATA_DIR, "grid10km_cdl_5yr_prop_intensity.csv")

# Output files
final_csv_path   <- file.path(PROCESSED_DATA_DIR, "grid10km_final_risk_attributes.csv")
final_gpkg_path  <- file.path(PROCESSED_DATA_DIR, "grid10km_final_risk.gpkg")

# --- Read inputs ---
g10_sf   <- st_read(gpkg_grid_path, quiet = TRUE)            # geometry + grid10_id
weather  <- read_csv(weather_csv_path, show_col_types = FALSE)
cdl_feat <- read_csv(cdl_csv_path,     show_col_types = FALSE)

# --- Light sanity checks ---
stopifnot("grid10_id" %in% names(g10_sf))
stopifnot("grid10_id" %in% names(weather))
stopifnot("grid10_id" %in% names(cdl_feat))

# If there are accidental duplicates by grid10_id, collapse them (keeps first non-NA)
dedupe_first_non_na <- function(df) {
  df %>%
    group_by(grid10_id) %>%
    summarise(dplyr::across(everything(), ~ .x[which.max(!is.na(.x))]), .groups = "drop")
}

weather  <- dedupe_first_non_na(weather)
cdl_feat <- dedupe_first_non_na(cdl_feat)

# --- Choose/rename columns to avoid clashes ---
# Keep weather metrics you actually need (adjust as needed)
weather_keep <- weather %>%
  select(
    grid10_id,
    # keep lon/lat from ONE source only; comment out if you prefer geometry centroids
    lon, lat,
    gdd, cum_gdd, n_days,
    # add other weather features if you computed them (examples):
    summer_avg_temp,
    summer_avg_percip,
    summer_heavy_rainfall_days,
    summer_hottest_temp,
    winter_coldest_temp,
    winter_extreme_cold_days,
    winter_warm_day_count
  ) %>%
  # guard against duplicated column names in join later
  rename(with_weather_lon = lon, with_weather_lat = lat)

cdl_keep <- cdl_feat %>%
  select(
    grid10_id, lon, lat,
    starts_with("crop_code_"),
    starts_with("prop1500_"),
    wei_prop, wei_intensity
  ) %>%
  rename(with_cdl_lon = lon, with_cdl_lat = lat)

# --- Join attributes ---
final_attr <- weather_keep %>%
  left_join(cdl_keep, by = "grid10_id")

# Optional: choose a single lon/lat to keep (CDL vs weather vs geometry centroids)
# Here we prefer CDL lon/lat if present, otherwise weather lon/lat.
final_attr <- final_attr %>%
  mutate(
    lon = dplyr::coalesce(with_cdl_lon, with_weather_lon),
    lat = dplyr::coalesce(with_cdl_lat, with_weather_lat)
  ) %>%
  select(-with_cdl_lon, -with_cdl_lat, -with_weather_lon, -with_weather_lat)

# --- Join back to geometry for mapping ---
# Keep only needed attributes to avoid bloat
final_sf <- g10_sf %>%
  left_join(final_attr, by = "grid10_id")

# --- Save outputs ---
# Attributes only
write_csv(final_attr, final_csv_path)
# With geometry (ready to plot in QGIS/ggplot)
st_write(final_sf, final_gpkg_path, delete_layer = TRUE, quiet = TRUE)

cat("Done.\n",
    "Rows in final table: ", nrow(final_attr), "\n",
    "Saved:\n - ", final_csv_path, "\n - ", final_gpkg_path, "\n", sep = "")

