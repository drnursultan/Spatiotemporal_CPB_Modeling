# R Codes

This section contains **R scripts** used for identifying crop types and calculating potato intensity with spatial analysis. The choice of R for these tasks is based on its efficient handling of raster data, supported by specialized packages that outperform Python for certain operations.
<br>

## Key Features

**Crop Type Identification** <br>
	• Crop types are identified using small buffer zones around each field, allowing precise spatial analysis.

**Potato Intensity Calculation** <br>
	• Potato intensity is computed using a **1500-meter buffer zone**, accounting for the spatial movement potential of Colorado Potato Beetles (CPB).

 # Why R?

R was chosen for these operations due to its robust geospatial libraries, which are highly optimized for raster data processing. Packages like terra, sf, and exactextractr enable efficient extraction and analysis of data, outperforming Python in similar use cases.

## Workflow

**R Packages Used** <br>
	•	terra: For raster data manipulation. <br>
	•	sf: For handling spatial vector data. <br>
	•	exactextractr: For extracting raster data within specified buffer zones. <br>

 
## How It Works
	1.	Buffer Zone Creation:
	•	A 1500-meter buffer zone is created around each field location.
	2.	Data Extraction:
	•	The exactextractr package extracts data for crop pixels (crop code 43 for potatoes) within the buffer zone.
	3.	Intensity Calculation:
	•	The function calculates the proportion of potato pixels relative to the total area within the buffer zone.
