## Model Overview
With this code, we are adding:
- land use land cover
	total sum of different land use types
	last observed land use type

## Dependencies
- Google Earth Engine
- Python API
- geemap

install.packages(c("curl", "httr", "openssl", "tidyverse", # required for tidyverse in Fedora
                   "terra", "rstan", "cmdstanr", "ggplot2"))


## Data Processing

### 1_satellite

	Intakes MAPBIOMAS raw fire data

	Outputs:
		Masks:


		Land use rasters:
		
-	Total number of fires
-	How many years since regrowth event was the last fire observed
	Intakes MAPBIOMAS raw land use/land cover data
-	Time since last observation of each land use type
-	Total years under each land use type
-	Last observed land use type before regrowth.





### 2_categorical
	Shapefiles from IBGE 2022: https://brasil.mapbiomas.org/en/tabela-de-camadas/

### 3_climate
	Intakes: 
Areas de proteção ambiental 2021 - https://centrodametropole.fflch.usp.br/
Areas de proteção ambiental 2023 - CNUC MMA


### 4_unifying_data


From Ma et al 2023
- fragmentation for 2020


## Modelling

### regrowth_functions.r


### run_model.r

	Intakes the rasters in model_ready_rasters.
	Makes sure all are in the same extent and resampled to the same resolution.
	Extracts values into a single dataframe, central_df
		Central_df is passed into the model
