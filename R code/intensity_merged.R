rm(list = ls())

library(progress)
library(tidyverse)
library(terra)
library(sf)
library(exactextractr)
library(stringr)
library(dplyr)


# Step 1: Read the input dataset
fields <- read.csv("merged_data/final_merged_data.csv")

# Step 2: Extract unique fields by latitude and longitude
unique_fields <- fields %>%
  select(lng, lat) %>%
  distinct()

# Assign unique field IDs
unique_fields <- unique_fields %>%
  mutate(field_id = row_number())

# Convert to an `sf` object
unique_fields_sf <- unique_fields %>%
  st_as_sf(coords = c("lng", "lat"), crs = 4326, remove = FALSE)

# Step 3: Calculate Majority Crop Code in Buffer Zone
get_majority_crop <- function(fields_sf, buffer_size = 10, cdl_years = 2008:2023) {
  lapply(cdl_years, function(yr) {
    # Load CDL raster for the given year
    r <- str_glue("../Project_X/webData/CDL_{yr}_55.tif")
    message("Processing Year: ", yr)
    
    # Buffer the points
    buffered_fields <- fields_sf %>%
      st_transform(5070) %>%
      st_buffer(buffer_size)
    
    # Extract majority crop data
    exact_extract(
      rast(r), buffered_fields,
      include_cols = "field_id",
      coverage_area = TRUE,
      summarize_df = TRUE,
      fun = function(df) {
        df %>%
          group_by(field_id, value) %>%
          summarise(
            pixel_count = n(),  # Count the number of pixels for each crop code
            .groups = "drop"
          ) %>%
          mutate(
            total_pixels = sum(pixel_count),  # Total pixels in the buffer zone
            buffer_radius = buffer_size
          ) %>%
          slice_max(pixel_count, n = 1) %>%  # Get majority crop by pixel count
          select(field_id, crop_code = value, pixel_count, total_pixels, buffer_radius)
      }
    ) %>%
      as_tibble() %>%
      mutate(cdl_year = yr)
  }) %>%
    bind_rows()
}

# Step 4: Process Data
cdl_years <- 2006:2023
buffer_size <- 50

majority_crop <- get_majority_crop(unique_fields_sf, buffer_size = buffer_size, cdl_years = cdl_years)

majority_crop <- majority_crop %>%
  group_by(field_id, cdl_year) %>%
  slice_max(pixel_count, n = 1, with_ties = FALSE) %>%  # Keep only one row per field per year
  ungroup()


# ------------- INTENSITY ---------------

get_potato_intensity <- function(fields_sf, buffer_size = 10, cdl_years = 2008:2023) {
  lapply(cdl_years, function(yr) {
    # Load CDL raster for the given year
    r <- str_glue("../Project_X/webData/CDL_{yr}_55.tif")
    message("Processing Year: ", yr)
    
    # Buffer the points
    buffered_fields <- fields_sf %>%
      st_transform(5070) %>%
      st_buffer(buffer_size)
    
    # Extract potato intensity
    exact_extract(
      rast(r), buffered_fields,
      include_cols = "field_id",
      coverage_area = TRUE,
      summarize_df = TRUE,
      fun = function(df) {
        df %>%
          group_by(field_id, value) %>%
          summarise(
            pixel_count = n(),  # Count the number of pixels for each crop code
            .groups = "drop"
          ) %>%
          mutate(
            potato_pixels = sum(pixel_count[value == 43], na.rm = TRUE),  # Pixels with crop code 43
            total_pixels = sum(pixel_count),  # Total pixels in the buffer zone
            buffer_radius = buffer_size
          ) %>%
          summarise(  # Summarise to keep potato intensity and other metrics
            field_id = first(field_id),
            potato_pixels = first(potato_pixels),
            total_pixels = first(total_pixels),
            intensity = potato_pixels / total_pixels,
            buffer_radius = first(buffer_radius)
          )
      }
    ) %>%
      as_tibble() %>%
      mutate(cdl_year = yr)
  }) %>%
    bind_rows()
}

cdl_years = 2010:2023
buffer_size <- 1500

# Using the new potato intensity function
potato_intensity <- get_potato_intensity(unique_fields_sf, buffer_size = buffer_size, cdl_years = cdl_years)

# Inspect the result
head(potato_intensity)



# ------ SAVING to CSV document ----------

potato_intensity <- potato_intensity %>%
  left_join(unique_fields %>% select(field_id, lng, lat), by = "field_id")

majority_crop <- majority_crop %>%
  left_join(unique_fields %>% select(field_id, lng, lat), by = "field_id")

write.csv(potato_intensity, "merged_data/3r_potato_intensity_by_field.csv", row.names = FALSE)
write.csv(majority_crop, "merged_data/2r_majority_crop_data.csv", row.names = FALSE)


# -------------- visualization part --------

plot_field_buffer <- function(fields_sf, f_id, year, buffer_size = 1500) {
  # Select the field based on field ID
  site <- fields_sf %>%
    filter(field_id == !!f_id) %>%  # Use the passed f_id parameter
    st_transform(5070)
  
  # Create buffer around the field
  buffer <- site %>% st_buffer(buffer_size)
  
  # Load CDL raster for the specified year
  cdl_raster <- rast(str_glue("../Project_X/webData/CDL_{year}_55.tif")) %>%
    crop(buffer)
  
  # Plot the raster and the buffer zone
  plot(cdl_raster, main = paste("Buffer Zone for Field", f_id, "in Year", year))
  plot(st_geometry(site), add = TRUE, col = "red", pch = 20)
  plot(st_geometry(buffer), add = TRUE, border = "blue", lwd = 2)
}

# Find the field and year with maximum potato intensity
max_potato_intensity <- potato_intensity %>%
  filter(intensity == max(intensity, na.rm = TRUE)) %>%
  slice(1)  # In case of ties, select the first row

field_with_max_intensity <- max_potato_intensity$field_id
year_with_max_intensity <- max_potato_intensity$cdl_year
max_intensity <- max_potato_intensity$intensity

cat("Max intensity is:", max_intensity, "by field ID:", field_with_max_intensity, "in year:", year_with_max_intensity, "\n")

plot_field_buffer(unique_fields_sf, 
                  f_id = field_with_max_intensity, 
                  year = year_with_max_intensity, 
                  buffer_size = 1500)




min_potato_intensity <- potato_intensity %>%
  filter(intensity == min(intensity, na.rm = TRUE)) %>%
  slice(1)  # In case of ties, select the first row


field_with_min_intensity <- min_potato_intensity$field_id
year_with_min_intensity <- min_potato_intensity$cdl_year
min_intensity <- min_potato_intensity$intensity


cat("Min intensity is:", min_intensity, "by field ID:", field_with_min_intensity, "in year:", year_with_min_intensity, "\n")

plot_field_buffer(unique_fields_sf, 
                  f_id = field_with_min_intensity, 
                  year = year_with_min_intensity, 
                  buffer_size = 1500)


# -------- MERGING with Original FIELDS dataset ------

# Add crop_type to fields
fields <- fields %>%
  left_join(majority_crop %>% select(lat, lng, cdl_year, crop_code), 
            by = c("lat" = "lat", "lng" = "lng", "year" = "cdl_year")) %>%
  rename(crop_type = crop_code)  # Rename crop_code to crop_type

# Add intensity to fields
fields <- fields %>%
  left_join(potato_intensity %>% select(lat, lng, cdl_year, intensity), 
            by = c("lat" = "lat", "lng" = "lng", "year" = "cdl_year"))


# Function to calculate potato proportion for last 5 years
calculate_potato_proportion <- function(lat, lng, year, majority_crop) {
  # Filter rows for the same location and last 5 years
  crop_history <- majority_crop %>%
    filter(lat == !!lat, lng == !!lng, cdl_year <= year, cdl_year > (year - 5))
  
  # Calculate the proportion of potato (crop_code == 43)
  potato_years <- sum(crop_history$crop_code == 43, na.rm = TRUE)
  total_years <- 5  # Fixed for 5-year window
  return(potato_years / total_years)
}

# Add potato_proportion to fields
fields <- fields %>%
  rowwise() %>%  # Apply the function row by row
  mutate(
    potato_proportion = calculate_potato_proportion(lat, lng, year, majority_crop)
  ) %>%
  ungroup()  # Remove row-wise grouping

# Save the updated fields dataset
write.csv(fields, "merged_data/4r_updated_final_fields.csv", row.names = FALSE)

# ------- BOX PLOT for intensity -------

library(ggplot2)

# Boxplot for potato intensity by year
potato_intensity %>%
  ggplot(aes(x = as.factor(cdl_year), y = intensity, group = cdl_year)) +
  geom_boxplot(fill = "lightblue") +
  labs(
    title = "Annual Potato Intensity",
    x = "Year",
    y = "Potato Intensity"
  ) +
  theme_minimal()

