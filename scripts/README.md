


## Model Overview


## Dependencies
- Google Earth Engine
- Python API
- geemap

## Scripts
### data_processing.py
	Shapefiles from IBGE 2022: https://brasil.mapbiomas.org/en/tabela-de-camadas/

Steps taken to reduce processing time

masked all the data to the age output before exporting



1_history_fire
	Intakes MAPBIOMAS raw fire data (./mapbiomas/lulc) and either location_forest_age.tif or santoro regrowth data ---- generates fire brick (location_fire_history.tif):
-	Total number of fires
-	How many years since regrowth event was the last fire observed

1_history_landuse
	Intakes MAPBIOMAS raw land use/land cover data (./mapbiomas/lulc) and 	location_forest_age.tif  ---- generates land use raster brick (location_lulc_history.tif):
-	Time since last observation of each land use type
-	Total years under each land use type
-	Last observed land use type before regrowth.

1_environmental
	Intakes: 
Areas de proteção ambiental 2021 - https://centrodametropole.fflch.usp.br/
Areas de proteção ambiental 2023 - CNUC MMA




### regrowth_functions.r


### run_model.r

	Intakes the rasters in model_ready_rasters.
	Makes sure all are in the same extent and resampled to the same resolution.
	Extracts values into a single dataframe, central_df
		Central_df is passed into the model


## References

MCWD data processing
Mapbiomas
