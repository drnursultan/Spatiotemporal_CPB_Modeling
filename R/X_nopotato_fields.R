# ============================================================
# File: test_nonpotato_finaldata_plots.R
# Purpose: Plot sample fields where croptype != 43 (non-potato)
# Uses: ../data/processed/final_data_for_modeling.csv
#       ../data/cdl/CDL_{year}_55.tif
# Outputs: ../fig/check_nonpotato_final/*.png
# ============================================================

library(tidyverse)
library(sf)
library(terra)
library(stringr)

# ----------------------------
# Paths (adjust if needed)
# ----------------------------
PROCESSED_DATA_DIR <- "../data/processed/"
CDL_DIR            <- "../data/cdl/"
OUT_DIR            <- "../fig/check_nonpotato_final/"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Read the final modeling data
# ----------------------------
fields <- read.csv(file.path(PROCESSED_DATA_DIR, "final_data_for_modeling.csv"))
message("fields: ", paste(dim(fields), collapse = " x "))

# Keep only the columns we need
df <- fields %>%
  select(farm, lat, lng, croptype, year)

# Quick sanity check on croptype codes
message("Unique croptype codes and counts:")
print(df %>% count(croptype, sort = TRUE))

# CDL potato class (typical): 43
POTATO_CODE <- 43

# Filter non-potato rows
non_potato <- df %>% filter(!is.na(croptype), croptype != POTATO_CODE)
message("Non-potato rows: ", nrow(non_potato))

# If nothing to plot, bail early
if (nrow(non_potato) == 0) {
  stop("No non-potato rows found (croptype != 43). Check your codes.")
}

# Make sf points (WGS84)
non_potato_sf <- non_potato %>%
  st_as_sf(coords = c("lng", "lat"), crs = 4326, remove = FALSE)

# ----------------------------
# Plot helper using coordinates directly
# ----------------------------
plot_point_buffer <- function(lon, lat, year, buffer_m = 300,
                              raster_dir = CDL_DIR,
                              main_title = NULL) {
  # 1) Create sf point in WGS84
  pt <- st_as_sf(data.frame(lon = lon, lat = lat),
                 coords = c("lon", "lat"), crs = 4326)
  
  # 2) Project to EPSG:5070 and buffer in meters
  pt_5070 <- st_transform(pt, 5070)
  buf_5070 <- st_buffer(pt_5070, buffer_m)
  
  # 3) Load CDL for that year
  rpath <- str_glue("{raster_dir}CDL_{year}_55.tif")
  if (!file.exists(rpath)) {
    warning("Missing raster for year ", year, ": ", rpath)
    return(invisible(NULL))
  }
  r <- rast(rpath)
  
  # 4) Transform buffer/point to raster CRS
  buf_rcrs <- st_transform(buf_5070, crs(r))
  pt_rcrs  <- st_transform(pt_5070, crs(r))
  
  # 5) Crop raster to buffer and plot
  r_crop <- crop(r, vect(buf_rcrs))
  
  if (is.null(main_title)) {
    main_title <- paste0("Non-potato | Year ", year)
  }
  
  plot(r_crop, main = main_title)
  plot(st_geometry(buf_rcrs), add = TRUE, border = "blue", lwd = 2)
  plot(st_geometry(pt_rcrs),  add = TRUE, col = "red", pch = 20, cex = 1.8)
}

# ----------------------------
# Pick first 10 non-potato rows to inspect
# ----------------------------
N <- min(10, nrow(non_potato_sf))
sample_rows <- non_potato_sf %>% slice_head(n = N)

# Loop: save PNGs
for (i in seq_len(nrow(sample_rows))) {
  lon  <- sample_rows$lng[i]
  lat  <- sample_rows$lat[i]
  yr   <- sample_rows$year[i]
  fm   <- sample_rows$farm[i]
  ct   <- sample_rows$croptype[i]
  
  fname <- file.path(
    OUT_DIR,
    paste0("nonpotato_", gsub("[^A-Za-z0-9]+","_", fm),
           "_y", yr, "_lon", round(lon,5), "_lat", round(lat,5), ".png")
  )
  
  message("Plotting: ", basename(fname), " (croptype=", ct, ")")
  
  png(fname, width = 1200, height = 1200, res = 150)
  plot_point_buffer(
    lon = lon, lat = lat, year = yr, buffer_m = 500, raster_dir = CDL_DIR,
    main_title = paste0(fm, " | Year ", yr, " | croptype=", ct)
  )
  dev.off()
}

message("Done. Check PNGs in: ", OUT_DIR)
