# Load required libraries
library(tidyverse)
library(terra)
library(sf)
library(exactextractr)
library(stringr)



# Step 1: Read the input dataset
fields <- read.csv("output_data/4_filtered_annual_data.csv")
dim(fields)  # Check dimensions

# Extract unique field IDs with their coordinates
unique_fields <- fields %>%
  select(field_fvid, lng, lat) %>%
  distinct()  # Retain only unique rows

# Convert to an `sf` object
unique_fields_sf <- unique_fields %>%
  st_as_sf(coords = c("lng", "lat"), crs = 4326, remove = FALSE)

print(dim(unique_fields_sf))  # Confirm the number of unique fields
head(unique_fields_sf)



# Step 3: Identify crop type with majority crop code from 49 pixels
# This one not works so good, I decided to calc CROP TYPE it through the Python ---------

get_majority_crop <- function(fields_sf, buffer_size = 10, cdl_years = 2010:2023) {
  lapply(cdl_years, function(yr) {
    # Load CDL raster for the given year
    r <- str_glue("../Project_X/webData/CDL_{yr}_55.tif")
    message("Processing Year: ", yr)
    
    # Buffer the points
    buffered_fields <- fields_sf %>%
      st_transform(5070) %>%
      st_buffer(buffer_size)
    
    # Extract data with pixel-based majority calculation
    exact_extract(
      rast(r), buffered_fields,
      include_cols = "field_fvid",
      coverage_area = TRUE,
      summarize_df = TRUE,
      fun = function(df) {
        df %>%
          group_by(field_fvid, value) %>%
          summarise(
            pixel_count = n(),  # Count the number of pixels for each crop code
            total_area = sum(coverage_area, na.rm = TRUE),  # Total buffer area
            .groups = "drop"
          ) %>%
          slice_max(pixel_count, n = 1) %>%  # Get majority crop by pixel count
          mutate(buffer_radius = buffer_size) %>%
          select(field_fvid, crop_code = value, pixel_count, total_area, buffer_radius)
      }
    ) %>%
      as_tibble() %>%
      mutate(cdl_year = yr)
  }) %>%
    bind_rows()
}

# Step 4: Extract potato area and total area within buffer zone
get_potato_extract <- function(pts, buffer_size = 1500, id_col = "field_fvid", cdl_years) {
  shp <- pts %>%
    st_transform(5070) %>% # Convert to meters
    st_buffer(buffer_size)
  
  lapply(cdl_years, function(yr) {
    r <- str_glue("../Project_X/webData/CDL_{yr}_55.tif")
    message("Processing: ", r)
    exact_extract(
      x = rast(r), y = shp,
      include_cols = id_col,
      coverage_area = TRUE,
      summarize_df = TRUE,
      fun = function(df) {
        df %>%
          rename(class = value) %>%
          mutate(total_area = sum(coverage_area), .by = all_of(id_col)) %>%
          summarise(
            class_area = sum(coverage_area[class == 43], na.rm = TRUE),
            total_area = total_area[1],
            intensity = class_area / total_area,
            .by = all_of(id_col)
          )
      }) %>%
      as_tibble() %>%
      mutate(cdl_year = yr, buffer_radius = buffer_size)
  }) %>%
    bind_rows()
}

# Step 5: Calculate 5-year intensity and potato proportion
calculate_potato_metrics <- function(potato_data, crop_codes, cdl_years) {
  potato_intensity_5yr <- lapply(cdl_years, function(yr) {
    potato_data %>%
      filter(cdl_year <= yr) %>%
      mutate(years_ago = yr - cdl_year) %>%
      filter(years_ago < 5) %>%
      group_by(field_fvid) %>%
      summarise(
        potato_5y_intensity = mean(intensity, na.rm = TRUE),
        potato_5y_weighted_intensity = weighted.mean(intensity, 5 - years_ago, na.rm = TRUE),
        potato_proportion = mean(crop_code == 43, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(year = yr)
  }) %>% bind_rows()
  return(potato_intensity_5yr)
}

# Define the CDL crop classes
cdl_classes <- tibble(
  class = c(43),
  crop_name = c("Potato")
)

# Define the years of interest
cdl_years <- 2010:2023

# Call the function with the unique fields dataset
potato_data <- get_potato_extract(unique_fields_sf, buffer_size = 1500, id_col = "field_fvid", cdl_years = cdl_years)

dim(potato_data)  # Dimensions of the resulting data
head(potato_data)

# Save the result for further use
write.csv(potato_data, "output_data/potato_intensity_by_field.csv", row.names = FALSE)

# Majority crop code calculation
crop_codes <- get_majority_crop(unique_fields_sf, buffer_size = 10, cdl_years = cdl_years)
dim(crop_codes)
crop_codes[(crop_codes$field_fvid == 803771)&(crop_codes$cdl_year==2018), ]


# Save crop_codes to a CSV file
write.csv(crop_codes, "output_data/crop_codes.csv", row.names = FALSE)



# Step 10: Function to visualize buffer zone for a specific field and year
plot_field_buffer <- function(fields_sf, field_id, year, buffer_size = 1500) {
  site <- fields_sf %>%
    filter(field_fvid == field_id) %>%
    st_transform(5070)
  
  buffer <- site %>% st_buffer(buffer_size)
  
  # Load CDL raster for the given year
  cdl_raster <- rast(str_glue("../Project_X/webData/CDL_{year}_55.tif")) %>%
    crop(buffer)
  
  # Plot the raster and buffer
  plot(cdl_raster, main = paste("Buffer Zone for Field", field_id, "in Year", year))
  plot(st_geometry(site), add = TRUE, col = "red", pch = 20)
  plot(st_geometry(buffer), add = TRUE, border = "blue", lwd = 2)
}

# Step 11: Visualize one field
set.seed(42)
random_field <- sample(fields$field_fvid, 1)
plot_field_buffer(unique_fields_sf, field_id = random_field, year = 2018)

plot_field_buffer(unique_fields_sf, field_id = 62996, year = 2010)

# Step 12: Boxplot for Potato Intensity
final_data %>%
  ggplot(aes(x = as.factor(year), y = intensity, group = year)) +
  geom_boxplot(fill = "lightblue") +
  labs(
    title = "Annual Potato Intensity",
    x = "Year",
    y = "Potato Intensity"
  ) +
  theme_minimal()