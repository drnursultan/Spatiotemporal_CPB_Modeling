# R/compute_grid_climate_2021_simple.R

# --- Packages ---
libs <- c("sf","terra","prism","dplyr","tibble","stringr","lubridate","readr")
invisible(lapply(libs, function(p){ if (!requireNamespace(p, quietly = TRUE)) install.packages(p); library(p, character.only = TRUE) }))
options(timeout = 1200)

# --- Inputs/paths ---
INTERIM_DATA_DIR <- "../data/interim"
GRID_GPKG <- file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.gpkg")
OUT_CSV   <- file.path(INTERIM_DATA_DIR, "grid10km_climate_2021.csv")
OUT_GPKG  <- file.path(INTERIM_DATA_DIR, "grid10km_enriched_climate_2021.gpkg")

stopifnot(file.exists(GRID_GPKG))
g10_sf  <- sf::st_read(GRID_GPKG, quiet = TRUE)
g10_deg <- st_transform(g10_sf, 4326)
cent    <- st_centroid(g10_deg)
coords  <- st_coordinates(cent)
ids     <- g10_sf$grid10_id
lon     <- coords[,1]; lat <- coords[,2]
pts4326 <- terra::vect(data.frame(lon=lon, lat=lat), geom=c("lon","lat"), crs="EPSG:4326")

# --- Time windows ---
YEAR <- 2021
gdd_start <- as.Date(sprintf("%d-04-01", YEAR))
gdd_end   <- as.Date(sprintf("%d-10-31", YEAR))
summer_start <- as.Date(sprintf("%d-05-01", YEAR))
summer_end   <- as.Date(sprintf("%d-08-31", YEAR))
winter_start <- as.Date(sprintf("%d-12-01", YEAR-1))
winter_end   <- as.Date(sprintf("%d-02-28", YEAR))
stopifnot(all(inherits(c(gdd_start,gdd_end,summer_start,summer_end,winter_start,winter_end), "Date")))
TBASE <- 10  # °C

# --- PRISM setup (persistent cache) ---
prism_dir <- file.path("..","data","prism_dl")
dir.create(prism_dir, recursive = TRUE, showWarnings = FALSE)
prism::prism_set_dl_dir(prism_dir)

retry_prism <- function(type, start, end, tries=5, sleep_sec=3){
  for(i in seq_len(tries)){
    ok <- try(prism::get_prism_dailys(type=type, minDate=start, maxDate=end, keepZip=FALSE), silent=TRUE)
    if(!inherits(ok, "try-error")) return(invisible(TRUE))
    Sys.sleep(sleep_sec*i)
  }
  stop("PRISM download failed for ", type, " ", start, "–", end)
}

# Download each window in one shot (much simpler, avoids month loops)
retry_prism("tmin", gdd_start,   gdd_end)
retry_prism("tmax", gdd_start,   gdd_end)

retry_prism("tmin", summer_start, summer_end)
retry_prism("tmax", summer_start, summer_end)
retry_prism("ppt",  summer_start, summer_end)

retry_prism("tmin", winter_start, winter_end)
retry_prism("tmax", winter_start, winter_end)
retry_prism("ppt",  winter_start, winter_end)

prism::prism_get_dl_dir()
head(prism::prism_archive_ls())

# --- FIX: list + stack helpers (no slashes in the pattern) ---
list_prism_files <- function(type, start, end){
  arch <- prism::prism_archive_ls()  # e.g. "PRISM_tmin_stable_4kmD2_20210501_bil"
  arch <- arch[grepl(paste0("^PRISM_", type, "_"), arch)]  # keep only that type
  bil  <- prism::pd_to_file(arch)  # -> full paths to the .bil files
  
  wanted <- format(seq.Date(start, end, by = "day"), "%Y%m%d")
  stamp  <- stringr::str_extract(basename(bil), "\\d{8}")
  bil[stamp %in% wanted]
}

stack_prism <- function(type, start, end){
  files <- list_prism_files(type, start, end)
  if (!length(files)) {
    stop("No PRISM ", type, " files found between ", start, " and ", end,
         ". Check prism_get_dl_dir() and that the window was downloaded.")
  }
  terra::rast(files)
}

# --- Now build the stacks from your existing cache ---
r_tmin_gdd <- stack_prism("tmin", gdd_start, gdd_end)
r_tmax_gdd <- stack_prism("tmax", gdd_start, gdd_end)

r_tmin_su  <- stack_prism("tmin", summer_start, summer_end)
r_tmax_su  <- stack_prism("tmax", summer_start, summer_end)
r_ppt_su   <- stack_prism("ppt",  summer_start, summer_end)

r_tmin_wi  <- stack_prism("tmin", winter_start, winter_end)
r_tmax_wi  <- stack_prism("tmax", winter_start, winter_end)
r_ppt_wi   <- stack_prism("ppt",  winter_start, winter_end)

# --- Extract at centroids ---
tmin_gdd <- terra::extract(r_tmin_gdd, pts4326)[,-1]
tmax_gdd <- terra::extract(r_tmax_gdd, pts4326)[,-1]
tmin_su  <- terra::extract(r_tmin_su,  pts4326)[,-1]
tmax_su  <- terra::extract(r_tmax_su,  pts4326)[,-1]
ppt_su   <- terra::extract(r_ppt_su,   pts4326)[,-1]
tmin_wi  <- terra::extract(r_tmin_wi,  pts4326)[,-1]
tmax_wi  <- terra::extract(r_tmax_wi,  pts4326)[,-1]
ppt_wi   <- terra::extract(r_ppt_wi,   pts4326)[,-1]

# --- Compute features ---
# GDD (Apr–Oct)
tmean_gdd <- (tmin_gdd + tmax_gdd)/2
gdd_day   <- pmax(tmean_gdd - TBASE, 0)
cum_gdd   <- rowSums(gdd_day, na.rm = TRUE)
gdd_avg   <- rowMeans(gdd_day, na.rm = TRUE)
n_days    <- rowSums(!is.na(gdd_day))

# Summer (May–Aug)
tmean_su <- (tmin_su + tmax_su)/2
summer_avg_temp            <- rowMeans(tmean_su, na.rm = TRUE)
summer_hottest_temp        <- apply(tmax_su, 1, max, na.rm = TRUE)
summer_avg_percip          <- rowMeans(ppt_su,  na.rm = TRUE)
summer_heavy_rainfall_days <- rowSums(ppt_su > 25, na.rm = TRUE)

# Winter (Dec prev – Feb)
tmean_wi <- (tmin_wi + tmax_wi)/2
winter_coldest_temp      <- apply(tmin_wi, 1, min, na.rm = TRUE)
winter_extreme_cold_days <- rowSums(tmin_wi <= -15, na.rm = TRUE)
winter_warm_day_count    <- rowSums(tmean_wi > 5, na.rm = TRUE)

# --- Assemble + save ---
clim_tbl <- tibble::tibble(
  grid10_id = ids,
  lon = lon, lat = lat,
  gdd = gdd_avg,
  cum_gdd = cum_gdd,
  n_days = n_days,
  summer_avg_temp            = summer_avg_temp,
  summer_avg_percip          = summer_avg_percip,
  summer_heavy_rainfall_days = summer_heavy_rainfall_days,
  summer_hottest_temp        = summer_hottest_temp,
  winter_coldest_temp        = winter_coldest_temp,
  winter_extreme_cold_days   = winter_extreme_cold_days,
  winter_warm_day_count      = winter_warm_day_count
)

readr::write_csv(clim_tbl, OUT_CSV)
message("Wrote: ", OUT_CSV)

g10_enriched <- g10_deg %>% left_join(clim_tbl, by = "grid10_id")
sf::st_write(g10_enriched, OUT_GPKG, delete_dsn = TRUE, quiet = TRUE)
message("Wrote: ", OUT_GPKG)

# --- Quick checks ---
print(head(clim_tbl))
print(summary(clim_tbl$cum_gdd, na.rm = TRUE))
