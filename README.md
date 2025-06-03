# Spatiotemporal Colorado Potato Beetle (CPB) Modeling

This project investigates the spatiotemporal dynamics of Colorado Potato Beetle (CPB) populations in Wisconsin. We integrate multi-source datasets to model CPB population trends and insecticide resistance risk.

👉 **Full project overview and process:**  
[View index.html](https://drnursultan.github.io/Spatiotemporal_CPB_Modeling/)

## Project Pipeline

1️⃣ **Data Extraction & Preprocessing (Python)**
- 01_extract_clean_combine.ipynb
- 02_add_climate_data.ipynb
- 03_outlier_multicollinearity_final_cleaning.ipynb

2️⃣ **Crop and Potato Features (R)**
- 01_crop_intensity_proportion.Rmd

3️⃣ **Final Dataset**
- data/processed/final_data_for_modeling.csv

4️⃣ **Modeling (R)**
- 02a_model_lmer_random_effects.Rmd
- 02b_model_gamm_spatial.Rmd
- 02c_model_spamm_spatial.Rmd

5️⃣ **EDA and Visualization (Python)**
- eda_visualization.ipynb

## Notes

- Large CDL raster files in `data/cdl/` are not uploaded due to GitHub size limits.
- Intermediate tables: `data/interim/`
- Final dataset: `data/processed/final_data_for_modeling.csv`

👉 **View the full report here:**  
🔗 [Open HTML Report](https://drnursultan.github.io/Spatiotemporal_CPB_Modeling/)
