# check_grid_centroids.R

# --- Packages ---
libs <- c("sf","dplyr","readr","ggplot2","tigris")
invisible(lapply(libs, function(p){
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))
options(tigris_use_cache = TRUE, tigris_class = "sf")

# --- Paths ---
INTERIM_DATA_DIR <- "../data/interim/"
csv_path  <- file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.csv")
gpkg_path <- file.path(INTERIM_DATA_DIR, "grid10km_enriched_noGDD.gpkg")

# --- 1) Read CSV with lon/lat to plot as red points ---
g10_attr <- read_csv(csv_path, show_col_types = FALSE)

# Make points from lon/lat (WGS84)
pts_wgs <- st_as_sf(g10_attr, coords = c("lon","lat"), crs = 4326)

# --- 2) Get Wisconsin boundary ---
wi_sf <- states(cb = TRUE) |>
  filter(STUSPS == "WI") |>
  st_transform(5070)   # Projected CRS for 10 km grid work

# --- 3) Get/Build 10 km grid polygons covering WI ---
if (file.exists(gpkg_path)) {
  message("Reading grid polygons from GPKG...")
  grid_sf <- st_read(gpkg_path, quiet = TRUE) |> st_transform(5070)
} else {
  message("GPKG not found; building a 10 km fishnet over WI…")
  # Build a 10 km fishnet, then clip to WI
  cell_km <- 10
  grid_sf <- st_make_grid(
    wi_sf,
    cellsize = c(cell_km*1000, cell_km*1000),
    square = TRUE,
    what = "polygons"
  ) |>
    st_as_sf() |>
    st_intersection(wi_sf) |>
    mutate(grid10_id = row_number())
}

# Ensure grid has an ID column for reference
if (!("grid10_id" %in% names(grid_sf))) {
  grid_sf <- grid_sf |> mutate(grid10_id = row_number())
}

# --- 4) Compute centroids of the grid polygons (to compare to CSV lon/lat) ---
grid_cent <- st_centroid(grid_sf) |>
  st_transform(4326) |>
  mutate(grid10_id = grid_sf$grid10_id)

# --- 5) Bring layers to a common CRS for plotting (WGS84) ---
wi_wgs   <- st_transform(wi_sf, 4326)
grid_wgs <- st_transform(grid_sf, 4326)

# --- 6) Plot: grid polygons + WI + CSV points + grid centroids ---
p <- ggplot() +
  geom_sf(data = wi_wgs, fill = "grey95", color = "white", linewidth = 0.3) +
  geom_sf(data = grid_wgs, fill = NA, color = "grey40", linewidth = 0.2) +
  # red points from CSV (what you want to verify)
  geom_sf(data = pts_wgs, color = "red", size = 1.6, alpha = 0.9) +
  # black X for centroids computed from polygons
  geom_sf(data = grid_cent, shape = 4, size = 1.6, stroke = 0.6, color = "black", alpha = 0.8) +
  coord_sf(expand = FALSE) +
  labs(title = "WI 10 km Grid: CSV lon/lat (red) vs. Polygon Centroids (black X)",
       subtitle = "If red points and black X’s overlap, CSV lon/lat are true centroids.",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

print(p)

# Optionally save
ggsave("../plots/maps/WI_grid_centroid_check.png", p, width = 9, height = 7, dpi = 300)
