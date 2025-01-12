# Spatiotemporal Colorado Potato Beetle (CPB) Modeling

This project investigates the spatiotemporal dynamics of Colorado Potato Beetle (CPB) populations in Wisconsin. By integrating multi-source datasets, we aim to uncover the factors influencing CPB distribution and abundance over time and space, supporting sustainable pest management practices.

## Project Overview

The dataset combines historical (2004–2018) and contemporary (2014–2023) CPB monitoring data from multiple research teams. Key data points include:
	•	Field-level information: Grower identities, field names, geolocations, and stage-specific CPB abundance metrics (adult, larval, and egg mass counts).
	•	Climate data: Sourced from MISTGEO, including precipitation, growing degree days (GDD), and temperature trends.
	•	Crop data: Extracted from the Cropland Data Layer (CDL) to calculate potato planting proportions and intensity.
	•	Potato metrics: Computed for each field over a 5-year span using R-based geospatial analysis.
 <br>

 ## Historical Dataset (2004–2018)
 
	•	Annual averages calculated for CPB abundance across all life stages.
	•	Climate data matched to field-level locations, with missing data imputed using nearest field data.
	•	Crop classifications and potato metrics integrated using spatial buffers.
	•	View Historical Dataset Details

 <br>

 ## Contemporary Dataset (2014–2023)
	•	Includes extended variables such as cumulative GDD, seasonal temperature extremes, and expanded field records.
	•	Computation methods align with the historical dataset for seamless integration.
	•	View Contemporary Dataset Details

 <br>

 ## Merged Dataset
	•	Historical and contemporary datasets merged with tolerance checks for field name changes and location shifts.
	•	Unmatched fields treated as new records for a comprehensive dataset.
	•	View Merged Dataset

 <br>

 ## Modeling and Analysis

The prepared dataset serves as the foundation for various predictive models and analytical insights:
	1.	Objectives:
	•	Predict insecticide resistance evolution using CPB resistance, genomic, abundance, and environmental data.
	•	Develop risk maps that integrate eco-physiological, weather, and potato intensity data.
	2.	Key Covariates:
	•	Climatic data: Monthly and annual precipitation and temperature.
	•	Potato intensity: Proportion and intensity metrics from CDL.
	•	CPB abundance: Annual weighted averages for different life stages.
	•	Insecticide resistance: Phenotypic measurements linked to imidacloprid resistance.
 ## Implemented Models
	•	Statistical and machine learning techniques were evaluated to find the best fit for CPB population predictions.
	•	Explore Models and Techniques
 
