

1_history_regrowth
Intakes MAPBIOMAS raw forest cover data (./mapbiomas/regrowth)------ generates forest_age raster:
-	Masks the history rasters with location_regrowth_mask.tif
-	Cleans the data:
o	Fixes typos
o	Standardizes redundant naming (all values within range 100-200 are anthropic, for example, and for “forest age” purposes all are designated the value ‘1’)
o	Removes unflagged regrowth moments
-	Returns raster with ALL SECONDARY FOREST AGES by using the identifier “503” as an indicator of regrowth (as indicated in MAPBIOMAS indexes) (location_forest_age.tif)

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
	WORLDCLIM tmin, tmax, and prec (2.5m resolution, stored in Polaris)
	FAO soil
	Output:
	Soil.rds
	prec_BRA_mean.tif
	prec_BRA_sd.tif
	temp_BRA_mean.tif
	temp_BRA_sd.tif

2_run_model

	Intakes the rasters in model_ready_rasters.
	Makes sure all are in the same extent and resampled to the same resolution.
	Extracts values into a single dataframe, central_df
		Central_df is passed into the model
