# --- Libraries ---
library(sf)
library(terra)
library(dplyr)
library(purrr)
library(stringr)
library(exactextractr)  # for polygon proportions (fast/robust)

# --- Paths ---
INTERIM_DATA_DIR  <- "../data/interim/"
CDL_DIR           <- "../data/cdl"   # files like CDL_2021_55.tif

# --- Load 10 km grid ---
g10 <- st_read(file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.gpkg"), quiet = TRUE)

# --- Centroids in WGS84 for lon/lat columns ---
cent_nat <- st_centroid(g10)          # native CRS (likely projected)
cent_4326 <- st_transform(cent_nat, 4326)

coords    <- st_coordinates(cent_4326)

g10_pts <- g10 %>%
  st_drop_geometry() %>%
  transmute(
    grid10_id,
    lon = coords[,1],
    lat = coords[,2]
  )
# attach centroid geometry (points, EPSG:4326)
st_geometry(g10_pts) <- st_geometry(cent_4326)
st_crs(g10_pts)      <- st_crs(cent_4326)

# --- 1.5 km buffers around centroids (in meters CRS) ---
cent_5070   <- st_transform(g10_pts, 5070)                 # CDL is in EPSG:5070 (US Albers)
buf_5070    <- st_buffer(st_geometry(cent_5070), 1500)     # 1.5 km
g10_buf_sf  <- st_sf(st_drop_geometry(cent_5070), geometry = buf_5070)
st_crs(g10_buf_sf) <- st_crs(cent_5070)                    # ensure CRS carried

# --- Years & weights ---
YEARS  <- 2021:2017
WTS    <- c(`2021`=5, `2020`=4, `2019`=3, `2018`=2, `2017`=1)
POTATO <- 43  # CDL crop code for potato

# --- Helper to load CDL raster for a year ---
load_cdl <- function(yr){
  f <- file.path(CDL_DIR, sprintf("CDL_%d_55.tif", yr))
  if (!file.exists(f)) stop("Missing CDL file: ", f)
  rast(f)  # terra SpatRaster, CRS ~ EPSG:5070 (NAD83 / Conus Albers)
}

# Base ID tables to enforce one row per grid
ids_pts  <- g10_pts      |> st_drop_geometry() |> dplyr::select(grid10_id)
ids_buf  <- g10_buf_sf   |> st_drop_geometry() |> dplyr::select(grid10_id)

# Make extractors always return (grid10_id, value) and then right_join to the IDs
get_centroid_crop <- function(r_cdl, pts_sf, buf_m = 50) {
  # Reproject points to raster CRS
  pts_rcrs  <- st_transform(pts_sf, crs(r_cdl))
  # Build small buffers
  small_buf <- st_buffer(st_geometry(pts_rcrs), buf_m)
  # Construct an sf that DEFINITELY has grid10_id as a column
  small_sf  <- st_sf(
    grid10_id = st_drop_geometry(pts_rcrs)$grid10_id,
    geometry  = small_buf,
    crs       = st_crs(pts_rcrs)
  )
  
  ee <- exactextractr::exact_extract(
    r_cdl, small_sf,
    include_cols = "grid10_id",
    summarize_df = TRUE,
    fun = function(df){
      if (nrow(df) == 0) return(tibble::tibble(crop_code = NA_integer_))
      df %>%
        dplyr::count(value, wt = coverage_fraction, name = "w") %>%
        dplyr::slice_max(w, n = 1, with_ties = FALSE) %>%
        dplyr::transmute(crop_code = as.integer(value))
    }
  ) %>% tibble::as_tibble()
  
  # If grid10_id didn’t come back, attach it explicitly (1 row per feature)
  if (!"grid10_id" %in% names(ee)) {
    ee <- dplyr::mutate(ee, grid10_id = small_sf$grid10_id, .before = 1)
  }
  
  # Sanity check
  stopifnot(nrow(ee) == nrow(pts_sf))
  dplyr::right_join(
    g10_pts %>% st_drop_geometry() %>% dplyr::select(grid10_id),
    ee,
    by = "grid10_id"
  )
}

get_buffer_potato_prop <- function(r_cdl, buf_sf){
  # Reproject buffers to raster CRS if needed
  if (st_crs(buf_sf) != st_crs(r_cdl)) buf_sf <- st_transform(buf_sf, crs(r_cdl))
  
  # Make sure grid10_id exists on the sf passed to exact_extract
  buf_sf <- buf_sf %>%
    dplyr::mutate(grid10_id = .data$grid10_id)  # no-op but ensures column is present
  
  res <- exactextractr::exact_extract(
    r_cdl, buf_sf,
    include_cols = "grid10_id",
    summarize_df = TRUE,
    fun = function(df){
      if (nrow(df) == 0) return(tibble::tibble(prop_potato = NA_real_))
      total <- sum(df$coverage_fraction, na.rm = TRUE)
      pot   <- sum(df$coverage_fraction[df$value == POTATO], na.rm = TRUE)
      tibble::tibble(prop_potato = ifelse(total > 0, pot/total, NA_real_))
    }
  ) %>% tibble::as_tibble()
  
  if (!"grid10_id" %in% names(res)) {
    res <- dplyr::mutate(res, grid10_id = buf_sf$grid10_id, .before = 1)
  }
  
  stopifnot(nrow(res) == nrow(buf_sf))
  dplyr::right_join(
    g10_buf_sf %>% st_drop_geometry() %>% dplyr::select(grid10_id),
    res,
    by = "grid10_id"
  )
}


centroid_list <- list()
buffer_list   <- list()

for (yr in YEARS){
  message("Processing CDL ", yr, " …")
  r <- load_cdl(yr)
  
  ctab <- get_centroid_crop(r, g10_pts, buf_m = 50) |>
    dplyr::rename(!!paste0("crop_code_", yr) := crop_code)
  
  btab <- get_buffer_potato_prop(r, g10_buf_sf) |>
    dplyr::rename(!!paste0("prop1500_", yr) := prop_potato)
  
  centroid_list[[as.character(yr)]] <- ctab   # has grid10_id + crop_code_yr
  buffer_list[[as.character(yr)]]   <- btab   # has grid10_id + prop1500_yr
}

# Base table
out <- g10_pts |>
  sf::st_drop_geometry() |>
  dplyr::select(grid10_id, lon, lat)

# Join yearly columns (now guaranteed to match ids)
out <- purrr::reduce(centroid_list, ~ dplyr::left_join(.x, .y, by = "grid10_id"), .init = out)
out <- purrr::reduce(buffer_list,   ~ dplyr::left_join(.x, .y, by = "grid10_id"), .init = out)




# --- Compute weighted 5-yr proportion & intensity ---
# centroid 0/1 potato indicator each year:
for (yr in YEARS) {
  code_col <- paste0("crop_code_", yr)
  bin_col  <- paste0("potato_", yr)
  out[[bin_col]] <- ifelse(out[[code_col]] == POTATO, 1, 0)
}

# Weighted proportion (centroid 0/1) : sum(w*y)/sum(w)
num_prop <- Reduce(`+`, lapply(YEARS, function(yr) WTS[as.character(yr)] * out[[paste0("potato_", yr)]]), init = 0)
den_prop <- sum(unname(WTS))
out$wei_prop <- num_prop / den_prop  # in [0,1]

# Weighted intensity (buffer proportion) : sum(w*prop1500)/sum(w)
num_int <- Reduce(`+`, lapply(YEARS, function(yr) WTS[as.character(yr)] * out[[paste0("prop1500_", yr)]]), init = 0)
den_int <- sum(unname(WTS))
out$wei_intensity <- num_int / den_int  # in [0,1]

# --- Final columns to keep ---
keep_cols <- c("grid10_id", "lon", "lat",
               paste0("crop_code_", YEARS),
               paste0("prop1500_", YEARS),
               "wei_prop", "wei_intensity")

out_final <- out %>% select(all_of(keep_cols))

# --- Save ---
write.csv(out_final, file.path(INTERIM_DATA_DIR, "grid10km_cdl_5yr_prop_intensity.csv"), row.names = FALSE)

# Quick peek
head(out_final)
