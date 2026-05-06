# Project Overview

We modelled the age of secondary forests in Brazil based on satellite data. This script imports and processes remote sensing data from Google Earth Engine, and compares different models to predict the biomass of secondary forests in 2020 and in the future.


```
forest_regrowth
├── 0_data
├── 1_gee
│   ├── 1_categorical.ipynb
│   ├── 2_edges_areas.ipynb
│   ├── 3_grids.ipynb
│   ├── 4_climate_soil.ipynb
│   ├── 5_land_use.ipynb
│   ├── 6_mature.ipynb
│   ├── 7_write_csv.ipynb
│   ├── 8_field_data.ipynb
│   ├── 9_projections.ipynb
│   ├── 10_extended_data.ipynb
│   └── utils.py
|
├── 2_modelling
│   ├── 0_groa_field_data.r
│   ├── 1_data_processing.r
│   ├── 1_parameters.r
│   ├── 2_cross_validate.r
│   ├── 2_forward_selection.r
│   └── 2_modelling.r
|
├── 3_analysis
│   ├── 0_asymptote_land_use.r
│   ├── 0_field.r
│   ├── 1_feature_importance.r
│   ├── 2_lag_field.r
│   ├── 3_future predictions.r
│   ├── 4_edge_biomass_hist.r
│   ├── 4_mature_distance_edge.r
│   ├── 4_pred_vs_obs_satellite
|   └── 4_tmf_comparison
│
├── README.md
└── requirements.txt
```


# 1_gee/

Scripts 1-6 process and export the data that that is then used to generate the dataframe for analysis in `7_write_csv`.



## 1_categorical.ipynb:
Exports images with binary masks for protected areas and indigenous land, or byte values for ecoregion and biome.

* **Imports:**
  * Indigenous land from [FUNAI](https://www.gov.br/funai/pt-br/atuacao/terras-indigenas/geoprocessamento-e-mapas)
  * Protected areas from [CEM - USP (2020)](https://centrodametropole.fflch.usp.br/pt-br/download-de-dados)
  * Biome data from [IBGE](https://www.ibge.gov.br/geociencias/informacoes-ambientais/vegetacao/15842-biomas.html)
  * Ecoregion from [RESOLVE](https://developers.google.com/earth-engine/datasets/catalog/RESOLVE_ECOREGIONS_2017)

* **Exports:**
  * `categorical` to GEE Image (indigenous, protected areas, biome)
  * `ecoreg` to GEE Image
  * `distance_to_border_mask` to GEE Image (pixels within 10km of a biome boundary—removing areas where the distance to nearest mature could be misinterpreted due to not including forests outside of Brazil, or forests of a different biome)



## 2_edges_areas.ipynb
Creates mask to consider for analysis only the secondary forest pixels that are surrounded by other secondary forest pixels on all sides (to avoid edge effects and biomass underestimations).

Estimates area of secondary forests per 1km². This is used to make total carbon sequestration estimates/projections.

* **Imports:**
  * Collection 9 MapBiomas Secondary Vegetation Age

* **Exports:**
  * `distance_to_secondary_edge` to GEE Image
  * `secondary_area_1km` to GEE Image
  * `pastureland_area_1km` to GEE Image



## 3_grids
To ensure proper spatial coverage while sparing compute time, we sample one pixel per 100km² grid cell to fit the model.

For final predictions, to ensure local specificity, we extract one pixel per 1km².

Also, for final carbon accumulation predictions, we select 10 points per cell in the Bezerra et al. 2022 dataset of future land use change conditions.

These grids are used in 7_write_csv to export the final dataframe for analysis.

* **Imports:**
  * Collection 9 MapBiomas Secondary Vegetation Age
  * Collection 9 MapBiomas Land Use Land Cover
  * `categorical`
  * `distance_to_border_mask`
  * `distance_to_secondary_edge`

* **Exports:**
  * `grid_10k_amazon_secondary_edge_removed` to GEE Feature Collection
  * `grid_10k_atlantic_secondary_edge_removed` to GEE Feature Collection
  * `grid_10k_amazon_secondary` to GEE Feature Collection
  * `grid_1k_amazon_secondary` to GEE Feature Collection
  * `grid_1k_amazon_pastureland` to GEE Feature Collection
  * `grid_Bezerra_10_points` to GEE Feature Collection: Samples 10 points per 100km² grid cell from the Bezerra et al. 2022 future predictions dataset.



##  4_climate_soil.ipynb:
Cleans and exports TerraClim and SoilGrids data for analysis.

* **TerraClim (yearly from 1958 to 2019):**
  * Summed (raw data is monthly totals):
    * Soil Moisture (mm)
    * Evapotranspiration (mm)
    * Precipitation (mm)
    * Climate Water Deficit (mm)

  * Averaged (raw data is monthly averages):
    * Temperature (C)
    * Vapour Pressure Deficit (kPa)
    * Palmer Drought Severity Index (PDSI)

  * Converted to kWh/m²/year (Total solar energy received per square meter over a year):
    * Solar Radiation (W/m^2)

* **SoilGrids (Averaged from 0-30cm depth and converted to the appropriate units):**
  * Bulk Density
  * Cation Exchange Capacity
  * Clay Content
  * Coarse fragments (> 2 mm)
  * Nitrogen
  * Organic Carbon Density
  * Soil Organic Carbon Stock
  * pH
  * Sand Content
  * Soil Organic Carbon

* **Exports:**
  * `terraclim_1958_2019` to GEE Image
  * `soilgrids` to GEE Image



## 5_land_use.ipynb:
Calculates:
  * Last observed land use type before regrowth
  * Sum of years under each land use type before regrowth
  * Number of fallow years

* **Imports:**
  * Collection 9 MapBiomas Land Use Land Cover

* **Exports:**
  * `land_use_non_aggregated_10yr`: Land use history restricted to 10 years from first to last observation of anthropogenic land use
  * `land_use_non_aggregated_5yr`: Land use history restricted to 5 years from first to last observation of anthropogenic land use
  * `land_use_non_aggregated_all`: Unrestricted land use history (All land use categories considered)
  * `land_use_aggregated_all`: Unrestricted land use history (Aggregated land use categories into 4 classes: pasture, perennial crops, annual crops, and mosaic)



## 6_mature.ipynb:
Process mature forest data to obtain surrounding mature forest cover and mature forest biomass estimates for the asymptote.

* **Imports:**
  * Collection 9 MapBiomas Land Use Land Cover
  * ESA CCI Biomass
  * MapBiomas forest age
  * Ecoregion
  * Amazon quarters (Heinrich et al 2021)

* **Exports:**
  * `distance_to_forest_edge`: distance to mature forest edge. Intermediate step to `distance_to_deep_forest`.
  * `sur_cover`: surrounding mature forest cover.
  * `distance_to_deep_forest`: distance to mature forests that are 1km or deeper.
  * `nearest_mature`: nearest mature forest biomass (10km x 10km grid averaged)
  * `quarters_ecoreg_biomass`: average biomass per Amazon quarter (from Heinrich et al 2021), used for Figure 2.



## 7_write_csv.ipynb:
Imports all data previously generated and exports it as CSV files for analysis.

Allows for the inclusion or exclusion of edge pixels (not surrounded by secondary forest on all sides) into the analysis.

Uses the projections from [Bezerra et al. 2022](https://doi.org/10.1371/journal.pone.0256052) to predict the future biomass accumulated by secondary forests.

* **Imports:**
  * Fire from MapBiomas Collection 3
  * Floodable Forests from MapBiomas Collection 9 (dummy variable)
  * Topography from ALOS
  * forest cover predictions for SSP1, SSP2 and SSP3 for 2015 and 2050.
  * `grid_10k_amazon_secondary_edge_removed` to GEE Feature Collection
  * `grid_10k_atlantic_secondary_edge_removed`
  * `grid_10k_amazon_secondary`
  * `grid_1k_amazon_secondary`
  * `grid_1k_amazon_pastureland`
  * `grid_Bezerra_10_points`

* **Exports:** CSV dataframes to run the model.
  * as a sample of 1 pixel per 100km² grid cell to ensure good spatial coverage without running into memory limits:
    * `grid_10k_amazon_secondary_edge_removed`
    * `grid_10k_amazon_secondary`: includes edges as a binary mask, to compare the biomass and conditions of the edge pixels
    * `grid_10k_atlantic_secondary_edge_removed`: samples the Atlantic Forest for a comparison of the land use predictive power outside of the Amazon

  * as a sample of 1 pixel per 1km² grid cell for future projections:
    * `grid_1k_amazon_pastureland`
    * `grid_1k_amazon_secondary`

  * `grid_Bezerra_10_points`: exports the regrowth conditions in areas expected to regrow by 2050 in at least one of the SSP scenarios. Exports the expected area to regrow per point.

  * `field_predictors.csv`: The predictors for the locations of the field data, which are then used in `3_analysis/0_field.r` to predict the biomass accumulation in the validation plots.



## 8_visualization.ipynb
Exports dataframes and maps that are not used for analysis, and are directly used to generate figures.

* **Figures:**
  * `Figure 2`: Asymptotes per aggregation level
  * `Figure 4`: Map of expected carbon accumulation by 2050
  * `Extended Data Fig. 1`: Biomass with EU TMF-obtained secondary forest age data
  * `Extended Data Fig. 2`: Mature Forest Biomass - distance to edge
  * `Extended Data Fig. 3`: Surrounding mature forest biomass across the Amazon
  * `Extended Data Fig. 4`: Why it is important to remove secondary forest edges



## utils.py:
  * defines project date range (1985-2020) and imports region of interest
  * defines export_image function
  * imports files for data processing (from Google Earth Engine)
  * select only pixels with exclusively the desired land use histories (exclude all instances of land use types we are not interested in)
  * makes grid cells for exporting data in gee_5_mature and gee_6_write_csv






# 2_modelling/

## 0_groa_field_data.r:
Imports and processes field data from the GROA project into a shapefile.
Shapefile is then used in 1_gee/8_field_data.ipynb for visualization
and to restrict the field data to the Amazon biome.

* **Imports:**
  * Field data from GROA (biomass_litter_CWD.csv)
  * site data from GROA (sites.csv)
* **Exports:**
  * Shapefile with aboveground biomass data for Brazil (field.shp)



## 1_data_processing.r:
Imports and processes the data for modelling.

* **Functions:**
  *`import_data`
    * Converts categorical variables to factors
    * Removes columns with extremely rare occurrences (less than 100 non-zero values)
    * Removes categorical values that occur less than 50 times
    * Selects which asymptote to use for the model (`nearest_mature`, `ecoreg_biomass`, `quarter_biomass`, `full_amazon`)
    * Splits the coordinates into a separate dataframe for export of results as shapefile
  *`normalize_independently`
    * Intakes training and testing dataframes
    * Normalizes the data for both dataframes based on the training data



## 1_parameters.r:
Defines the categories of parameters:
  * Land Use (lu, fallow, num_fires)
  * Soil
  * Categorical
  * Binary
  * Landscape (dist, sur_cover)
Defines the parameter lists for comparisons in `3_analysis.r`
  * basic_pars_options
  * data_pars_options



## 2_cross_validate.r:
  Evaluates the model performance using 5-fold cross-validation.
* **Functions:**
  *`calc_r2`
  *`cross_validate`



## 2_forward_selection.r:
  Iteratively fits parameter combinations with `run_optim()` and selects the one that minimizes the Akaike Information Criterion (AIC), excluding parameters that do not improve the model.

* **Functions:**
  *`find_combination_pars`



## 2_modelling.r:
Defines the main functions for the modelling process.

* **Functions:**
  *`run_optim`
  *`calc_rss`
  *`growth_curve`






# 3_analysis

## 0_asymptote_land_use.r:
Compares the R2 values of different models trained on satellite data.

Comparisons:
  * Asymptotes (`nearest_mature`, `ecoreg_biomass`, `quarter_biomass`, `full_amazon`)
  * Land Use (`non_aggregated_all`, `aggregated_all`, `non_aggregated_5yr`, `non_aggregated_10yr`)
  * Biomes (Amazon, Atlantic Forest)

* **Imports:** (all CSVs in directories)
  * `grid_10k_amazon_secondary`
  * `land_use_aggregated_all`
  * `land_use_non_aggregated_all`
  * `land_use_non_aggregated_5yr`
  * `land_use_non_aggregated_10yr`

* **Exports:**
  * `0_asymptotes.csv`: R2 values for each asymptote
  * `0_land_use.csv`: R2 values for each land use type per biome



## 0_field.r:
Obtains the R2 values for the field data based on the model trained from satellite data

* **Imports:**
  * `grid_10k_amazon_secondary`
  * `field_predictors.csv`
* **Exports:**
  * `0_field_results.csv`: R2 value for the field data and fit theta value from field data
  * `field_age_histogram.png`: Histogram of field data ages
  * `predicted_vs_observed_field.png`: Scatterplot of predicted vs observed biomass for field data




## 1_feature_importance.r:
Figure 2: Barplots.
Compares the relative importance of the parameters of full_amazon (inflexible) asymptote with the R2 of the nearest_mature (flexible) asymptote.
Compares R2 with three levels of asymptote aggregation with just age as the only predictor.
Compares the relative importance of the parameters of the Amazon and Atlantic Forest models (NN asymptote) to show land use is not incorporated.

* **Imports:** (all CSVs in directories)
  * `grid_10k_amazon_secondary`
  * `land_use_aggregated_all`
  * `land_use_non_aggregated_all`
  * `land_use_non_aggregated_5yr`
  * `land_use_non_aggregated_10yr`
* **Exports:**
  * figure_2_model_performance.jpg
  * figure_2_asymptote_barplot.jpg



## 2_lag_field.r:
Growth curve line graph.
Compares the growth rate of intercept and lag models.
Overlays the average biomass per age from the field data scatterplot.

* **Imports:**
  * `grid_10k_amazon_secondary`
  * `field_predictors.csv`
* **Exports:**
  * `lag_field_biomass.jpeg`
  * `lag_field_biomass_legend.jpeg`



## 3_future_predictions.r:
Barplot 1: Compares the biomass gain by 2050 for:
  * random 5% of pastureland
  * 5% with top regrowth potential
  * all secondary forests
Barplot 2: Shows current area of:
  * 5% of pastureland
  * secondary forests
Barplot 3: Shows current biomass stocked in:
  * random 5% of pastureland
  * secondary forests
Shapefile 1: Predicted biomass gain by 2050 for all pastureland.
Shapefile 2: Predicted biomass gain by 2050 for all secondary forests.

* **Imports:**
  * `grid_1k_amazon_secondary`: all CSVs in directory
  * `grid_1k_amazon_pastureland`: all CSVs in directory
* **Exports:**
  * `figure_4_c.jpeg`
  * `figure_4_d.jpeg`
  * `pred_2050_pastureland_all.shp`
  * `pred_2050_secondary_all.shp`



## 4_mature_distance_edge.r:
Shows the biomass of mature forests in relation to the distance to the nearest forest edge.

* **Imports:**
  * `mature_biomass_distance.csv`
* **Exports:**
  * `mature_biomass_distance_edge.jpeg`



## 4_pred_vs_obs_satellite.r:
Shows the predicted vs observed biomass for the satellite data.

* **Imports:**
  * `grid_10k_amazon_secondary`
* **Exports:**
  * `predicted_vs_observed_satellite.png`



## Getting Started
To get started with the project, follow these steps:

1. Clone the repository:
  ```sh
  git clone <repository-url>
  cd forest_regrowth
  ```

2. Install the required dependencies:
  ```sh
  pip install -r requirements.txt
  ```

Refer to `requirements.txt` for the complete list of dependencies.