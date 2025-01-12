# Spatiotemporal_CPB_Modeling
This project investigates the spatiotemporal dynamics of Colorado Beetle (CPB) populations. 


It explores data analysis and modeling techniques to understand the factors influencing CPB distribution and abundance over space and time.

Here’s a refined and detailed description for your README.md file, incorporating your project’s key aspects and using a polished style.

Spatiotemporal Colorado Potato Beetle (CPB) Modeling

This project investigates the spatiotemporal dynamics of Colorado Potato Beetle (CPB) populations in Wisconsin, leveraging an integrated dataset to understand the factors influencing their distribution and abundance over time and space. The study aims to support sustainable pest management practices through advanced data analysis and predictive modeling.

Project Overview

The dataset combines historical (2004–2018) and contemporary (2014–2023) CPB monitoring data from multiple research teams. It includes detailed field-level information such as grower identities, field names, precise geolocations, and stage-specific CPB abundance metrics (adult, larval, and egg mass counts). Additional variables are integrated to enrich the dataset:
	1.	Climate Data: Collected from MISTGEO, covering key variables like precipitation, growing degree days (GDD), and temperature extremes.
	2.	Crop Data: Extracted from the Cropland Data Layer (CDL) using USDA tools, providing annual crop classifications, potato planting proportions, and spatial potato intensity.
	3.	Potato Metrics: Calculated proportion and intensity metrics over a 5-year span for each field using R-based geospatial analysis.

Data Processing Workflow

Historical Dataset
	•	Covers the years 2004–2018.
	•	Annual averages were computed for CPB abundance across all life stages.
	•	Climate data was matched with field-level data based on year and location. Missing climate records were imputed using nearest available data.
	•	Crop type and potato planting data were extracted, including buffer-based calculations for fields within a 1.5 km radius.
	•	Potato proportion and intensity metrics were added, emphasizing multi-year planting trends.

Contemporary Dataset
	•	Extends data to 2023 with additional variables such as cumulative GDD and seasonal temperature trends.
	•	Data was processed with similar methods as the historical set, ensuring compatibility and consistency across datasets.

Merging Datasets
	•	Historical and contemporary datasets were merged, accommodating shifts in field names and slight relocations with spatial and temporal tolerance checks.
	•	Unmatched fields were treated as new records, ensuring a comprehensive dataset.

Modeling and Analysis

The final dataset was used to explore modeling approaches tailored to CPB population dynamics:
	•	Objectives:
	•	Predict insecticide resistance evolution using CPB resistance, genomic, abundance, and environmental data.
	•	Develop risk maps integrating eco-physiological insights, weather variability, and potato planting intensity.
	•	Key Covariates:
	•	Climatic variables: monthly and annual precipitation and temperature.
	•	Potato intensity: proportion and intensity metrics derived from CDL data.
	•	CPB abundance: annual weighted averages for multiple life stages.
	•	Insecticide resistance: measured phenotypes linked to imidacloprid resistance.

References

This project builds on methodologies and insights from leading studies in agricultural pest dynamics, including landscape ecology and risk mapping approaches ￼ ￼ ￼ ￼ ￼ ￼.

Would you like further assistance in linking specific files, adding citations, or including any additional sections?
