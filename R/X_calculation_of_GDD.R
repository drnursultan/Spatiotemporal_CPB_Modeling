# --- Libraries ---
library(sf)
library(terra)
library(prism)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
options(timeout = 1200)

# --- 0) Paths & grid centroids ---
INTERIM_DATA_DIR <- "../data/interim/"
g10_path <- file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.gpkg")

g10_sf <- st_read(g10_path, quiet = TRUE)
cent   <- st_transform(st_centroid(g10_sf), 4326)
cc     <- st_coordinates(cent)

ids <- g10_sf$grid10_id
lon <- cc[,1]
lat <- cc[,2]

# --- 1) Period & base ---
YEAR  <- 2021
START <- as.Date(sprintf("%d-04-01", YEAR))
END   <- as.Date(sprintf("%d-10-31", YEAR))
TBASE <- 10  # °C

# --- 2) Make sure PRISM download dir is right (you already have files there) ---
# This should print the temp dir you showed earlier:
print(prism_get_dl_dir())

# If you ever need to point to a permanent folder, use:
# prism_set_dl_dir("D:/PRISM")  # example

# (Re)download if needed (safe to skip if already present)
# Monthly chunks are more robust; comment out if you already downloaded.
dl_month <- function(type, y, m, tries = 3){
  from <- as.Date(sprintf("%04d-%02d-01", y, m))
  to   <- from + months(1) - days(1)
  for (k in 1:tries){
    ok <- try(get_prism_dailys(type = type, minDate = from, maxDate = to, keepZip = FALSE),
              silent = TRUE)
    if (!inherits(ok, "try-error")) return(invisible(TRUE))
    Sys.sleep(2 * k)
  }
  stop("PRISM download failed for ", type, " ", y, "-", m)
}
for (m in 4:10){
  dl_month("tmin", YEAR, m)
  dl_month("tmax", YEAR, m)
}

# --- 3) Find the exact daily .bil files in the archive for our date window ---
arch_dirs <- prism_archive_ls()  # character vector like "PRISM_tmin_stable_4kmD2_20210401_bil"

# helper to filter archive -> ordered file paths within date window
collect_daily_files <- function(type, start, end){
  # keep only the right product family
  keep <- arch_dirs[grepl(paste0("^PRISM_", type, "_"), basename(arch_dirs))]
  # convert archive dirs -> full .bil file paths
  files <- pd_to_file(keep)
  # extract YYYYMMDD stamp from filename
  stamp <- str_extract(basename(files), "\\d{8}")
  # keep files within [start, end]
  wanted <- format(seq.Date(start, end, by = "day"), "%Y%m%d")
  files  <- files[!is.na(stamp) & stamp %in% wanted]
  stamp  <- stamp[!is.na(stamp) & stamp %in% wanted]
  # order by date to align layers
  ord <- order(stamp)
  files[ord]
}

files_tmin <- collect_daily_files("tmin", START, END)
files_tmax <- collect_daily_files("tmax", START, END)

stopifnot(length(files_tmin) > 0, length(files_tmax) > 0,
          length(files_tmin) == length(files_tmax))

# --- 4) Stack rasters in matching order ---
r_tmin <- rast(files_tmin)
r_tmax <- rast(files_tmax)

# --- 5) Extract temps at centroids ---
pts <- vect(data.frame(lon = lon, lat = lat), geom = c("lon","lat"), crs = "EPSG:4326")

tmin_mat <- as.matrix(extract(r_tmin, pts)[,-1, drop = FALSE])
tmax_mat <- as.matrix(extract(r_tmax, pts)[,-1, drop = FALSE])

# --- 6) Compute daily mean and GDD, then seasonal summaries ---
tmean_mat <- (tmin_mat + tmax_mat)/2
gdd_day   <- pmax(tmean_mat - TBASE, 0)        # Type-B method (no upper cap)

cum_gdd <- rowSums(gdd_day, na.rm = TRUE)      # sum across days
n_days  <- rowSums(!is.na(gdd_day))            # count of valid days
gdd_avg <- cum_gdd / n_days                    # optional daily average

gdd_tbl <- tibble(
  grid10_id = ids,
  lon = lon, lat = lat,
  gdd = gdd_avg,
  cum_gdd = cum_gdd,
  n_days = n_days
)

# --- 7) (Optional) save & quick sanity check ---
# write.csv(gdd_tbl, file.path(INTERIM_DATA_DIR, "grid10km_gdd_2021_base10C.csv"), row.names = FALSE)

print(head(gdd_tbl))
summary(gdd_tbl$cum_gdd)

# keep only the new metrics for the join
gdd_keep <- gdd_tbl %>%
  select(grid10_id, gdd, cum_gdd, n_days)

# transform grid, compute centroids once, and join by ID
g10_out <- g10_sf %>%
  st_transform(4326) %>%
  mutate(
    lon = st_coordinates(st_centroid(st_geometry(.)))[,1],
    lat = st_coordinates(st_centroid(st_geometry(.)))[,2]
  ) %>%
  left_join(gdd_keep, by = "grid10_id")

# quick check
g10_out %>% st_drop_geometry() %>% select(grid10_id, lon, lat, gdd, cum_gdd, n_days) %>% head()

g10_out %>% 
  st_drop_geometry() %>% 
  write.csv("../data/interim/grid10km_with_gdd.csv", row.names = FALSE)

st_write(
  g10_out, 
  "../data/interim/grid10km_with_gdd.gpkg", 
  delete_dsn = TRUE
)
