# Spatiotemporal Colorado Potato Beetle (CPB) Modeling

This project investigates the spatiotemporal dynamics of Colorado Potato Beetle (CPB) populations in Wisconsin. We integrate long-term scouting data, climate records, and crop history information to model CPB abundance and insecticide resistance risk. The goal is to understand key ecological drivers and to provide predictive tools that can support regional pest management.

👉 **Full project overview and process:**  
[View index.html](https://drnursultan.github.io/Spatiotemporal_CPB_Modeling/)

## Project Pipeline

1️⃣ **Data Extraction & Preprocessing (Python)**  
Scripts for cleaning and preparing the scouting and weather data:
- `01_extract_clean_combine.ipynb`  
- `02_add_climate_data.ipynb`  
- `03_outlier_multicollinearity_final_cleaning.ipynb`

2️⃣ **Crop and Potato Features (R)**  
Calculation of potato proportion and intensity using CDL data:
- `01_crop_intensity_proportion.Rmd`

3️⃣ **Final Dataset**  
Processed dataset ready for modeling:
- `data/processed/final_data_for_modeling.csv`

4️⃣ **Modeling (R)**  
Comparisons of different statistical frameworks:
- `02a_model_lmer_random_effects.Rmd`  
- `02b_model_gamm_spatial.Rmd`  
- `02c_model_spamm_spatial.Rmd`  
- `risk_map_model.html` (reduced spaMM model used for statewide risk mapping)

5️⃣ **EDA and Visualization (Python)**  
Exploratory analysis and plotting:
- `eda_visualization.ipynb`

## Risk Maps

We extended the modeling framework to statewide **predictive risk mapping**.  
- A 10 × 10 km grid was created across Wisconsin, and centroids were used to assign predictors.  
- Climate predictors were derived from PRISM data, and potato intensity/proportion metrics were aggregated within each grid.  
- The reduced spaMM model (significant predictors + spatial correlation + year effects, excluding farm) was chosen for risk mapping due to stability and similar performance to the full model.  

Outputs include:
- Baseline 2021 risk surface  
- Scenario maps with reduced (×0.7) and increased (×1.5) potato intensity  

These maps illustrate how crop history and climate shape CPB risk across the state and serve as tools for proactive pest management.

## Notes

- Large CDL raster files in `data/cdl/` are not uploaded due to GitHub size limits.  
- Intermediate tables are stored in `data/interim/`.  
- Final modeling dataset is available in `data/processed/final_data_for_modeling.csv`.  
- Processed statewide risk surfaces are in `data/processed/` as `.csv` and `.gpkg`.

👉 **View the full report here:**  
🔗 [Open HTML Report](https://drnursultan.github.io/Spatiotemporal_CPB_Modeling/)