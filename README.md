## Project Overview
ReSELVA (Remote sensing SEcondary forest Local Variability Analysis) is a model that predicts the age of secondary forests in Brazil based on satellite data. This script imports and processes remote sensing data from Google Earth Engine, and compares different models to predict the biomass of secondary forests in 2020 and in the future.


```
forest_regrowth
├── 0_data
├── 1_gee
│   ├── 0_utils.py
│   ├── 1_categorical.ipynb
│   ├── 2_age_biomass.ipynb
│   ├── 3_grids_area.ipynb
│   ├── 4_climate_soil.ipynb
│   ├── 5_land_use.ipynb
│   ├── 6_mature.ipynb
│   ├── 7_write_csv.ipynb
│   └── 8_field_data.ipynb
|
├── 2_modelling
│   ├── 0_multicollinearity.r
│   ├── 0_groa_field_data.r
│   ├── 1_data_processing.r
│   ├── 1_parameters.r
│   ├── 2_cross_validate.r
│   ├── 2_feature_selection.r
│   └── 2_modelling.r
|
├── 3_analysis
│   ├── 0_asymptote_land_use.r
│   ├── 0_field.r
│   ├── 1_feature_importance.r
│   ├── 1_model_performance.r
│   ├── 2_lag_field.r
│   ├── 3_future predictions.r
|   └── EXT_pred_vs_obs.r
│
├── README.md
└── requirements.txt
```


### 1_gee_scripts/
- **gee_0_utils.py**:
    - defines project date range (1985-2020) and imports region of interest
    - defines export_image function
    - imports files for data processing (from Google Earth Engine)
    - select only pixels with exclusively the desired land use histories (exclude all instances of land use types we are not interested in)
    - makes grid cells for exporting data in gee_5_mature and gee_6_write_csv

- **gee_2_categorical.ipynb**:
    Imports:
      - Indigenous Lands
      - Protected Areas
      - Ecoregions
      - Biomes
    Exports:
      - "categorical" to GEE assets
      - "distance_to_border_mask" (pixels within 10km of a biome boundary - removing areas where the distance to nearest mature could be misinterpreted due to not including forests outside of Brazil, or forests of a different biome)

- **gee_2_age_biomass.ipynb**:
    - exports secondary forest age data from TMF and Mapbiomas
      - removes pixels with ages that don't match the IPCC estimates
      - removes isolated pixels (keeps only pixels within a patch of at least 1 hectare)
      - removes pixels within 10km of a biome boundary (distance_to_border_mask)

- **gee_3_climate_soil.ipynb**:
    - TerraClim
    Calculated yearly metrics.
        Summed:
        - Solar Radiation
        - Soil Moisture
        - Precipitation
        Averaged:
        - Temperature
        - Vapour Pressure
        - Evapotranspiration
    "yearly_terraclim": values of the metrics across all years from 1985-2019, and the means across time
    
    - SoilGrids
        - Bulk Density
        - Cation Exchange Capacity
        - Clay Content
        - Coarse fragments (> 2 mm)
        - Nitrogen
        - Organic Carbon Density
        - Soil Organic Carbon Stock
        - pH
        - Sand Content
        - Soil Organic Carbon
    All averaged from 0-30cm depth and converted to the correct units.
    
    - CMIP6 Climate
    Gets yearly climate data from the CMIP6 dataset for the period 1985-2019, and averages it across all years.
        - Near Surface Specific Humidity
        - Near Surface Air Temperature
        - Moisture in Upper Portion of Soil Column
        - Precipitation
        - Surface Downwelling Shortwave Radiation
    Exports:
        - "CMIP6_historical"
        - "CMIP6_ssp245"
        - "CMIP6_means"

- **gee_4_mature.ipynb**:
    Imports:
        - ESA CCI Biomass
        - MapBiomas forest age
        - TMF forest age
        - Ecoregions
        - Amazon quarters (Heinrich et al 2021)
    Exports:
        - "distance_to_forest_edge"
        - "sur_cover"
        - "mature_biomass"
        - "mature¨biomass_exclude_edge"
        - "nearest_mature_biomass"
        - "quarters_ecoreg_biomass"

- **gee_5_land_use.ipynb**:
    Imports land use Collection 9 data from MapBiomas.
    Calculates:
        - Last observed land use type before regrowth
        - Sum of years under each land use type before regrowth
        - Number of fallow years
    Five of such dataframes are exported:
        - Land use history restricted to 15 years from first to last observation of anthropogenic land use
        - Land use history restricted to 10 years from first to last observation of anthropogenic land use
        - Land use history restricted to 5 years from first to last observation of anthropogenic land use
        - Unrestricted land use history (All land use categories considered)
        - Unrestricted land use history (Aggregated land use categories into 4 classes: pasture, perennial crops, annual crops, and mosaic)

- **gee_6_write_csv.ipynb**:
    Imports all data previously generated and exports it as CSV files for analysis.
    - Biomass and age data source comparisons
    - Mature forest biomass comparisons
    - Field Data
    - Main model data assemblage
      - New introductions:
        - Fire from MapBiomas
        - Floodable Forests from MapBiomas (dummy variable)
        - Topography from ALOS
        Exports all data required for the model:
        - as a random sample of 15000 pixels per biome (for model fitting and validation)
        - as tiles incorporating all pixels classified as secondary forests for 2020 (for future predictions)


- **gee_8_field_data.ipynb**:
    Imports field.shp from 2_modelling/groa_field_data.r



### 2_modelling/

- **0_groa_field_data.r**:
    Imports and processes field data from the GROA project into a shapefile.
    Shapefile is then used in 1_gee/8_field_data.ipynb for visualization
    and to restrict the field data to the Amazon biome.
    Imports:
    - Field data from GROA (biomass_litter_CWD.csv)
    - site data from GROA (sites.csv)
    Exports:
    - Shapefile with aboveground biomass data for Brazil (field.shp)

- **0_multicollinearity.r**:
    Tests for multicollinearity with VIF. Multicollinear variables are removed in 1_parameters.r

- **1_data_processing.r**:
    Imports and processes the data for modelling.
    Functions:
    - "import_data"
      - Converts categorical variables to factors
      - Removes columns with extremely rare occurrences (less than 100 non-zero values)
      - Removes categorical values that occur less than 50 times
      - Selects which asymptote to use for the model ("nearest_mature", "ecoreg_biomass", "quarter_biomass", "full_amazon")
      - Splits the coordinates into a separate dataframe for export of results as shapefile
    - "normalize_independently"
      - Intakes training and testing dataframes
      - Normalizes the data for both dataframes based on the training data

- **1_parameters.r**:
    Defines the categories of parameters:
    - Land Use (lu, fallow, num_fires)
    - Soil
    - Categorical
    - Binary
    - Landscape (dist, sur_cover)
    Defines the parameter lists for comparisons in 3_analysis.r
    - basic_pars_options
    - data_pars_options

- **2_cross_validate.r**:
    Evaluates the model performance using 5-fold cross-validation.
    Functions:
    - "calc_r2"
    - "cross_validate"

- **2_forward_selection.r**:
    Iteratively fits parameter combinations with `run_optim()` and selects the one that minimizes the Akaike Information Criterion (AIC), excluding parameters that do not improve the model.
    Functions:
      - "find_combination_pars"

- **2_modelling.r**:
    Defines the main functions for modelling:
    - "run_optim"
    - "calc_rss"
    - "growth_curve"


### 3_analysis

- **0_asymptote_land_use.r**:
    Compares the R2 values of different models trained on satellite data.
    Comparisons:
    - Asymptotes ("nearest_mature", "ecoreg_biomass", "quarter_biomass", "full_amazon")
    - Land Use (non_aggregated_all, aggregated_all, non_aggregated_5yr, non_aggregated_10yr)
    - Biomes (Amazon, Atlantic Forest)

    Inputs:
    - "grid_10k_amazon_secondary": all CSVs in directory
    
    Outputs:
    - "0_asymptotes.csv": R2 values for each asymptote
    - "0_land_use.csv": R2 values for each land use type per biome

- **0_field.r**:
    Obtains the R2 values for the field data based on the model trained from satellite data
    Inputs:
    - grid_10k_amazon_secondary
    - field_predictors.csv
    Outputs:
    - "0_field_results.csv": R2 value for the field data and fit theta value from field data
    - "field_age_histogram.png": Histogram of field data ages
    - "predicted_vs_observed_field.png": Scatterplot of predicted vs observed biomass for field data

- **1_feature_importance.r**:
    Figure 2: Barplots.
    Compares the relative importance of the parameters of full_amazon (inflexible) asymptote with the R2 of the nearest_mature (flexible) asymptote.
    Compares R2 with three levels of asymptote aggregation with just age as the only predictor.
    Compares the relative importance of the parameters of the Amazon and Atlantic Forest models (NN asymptote) to show land use is not incorporated.

- **2_lag_field.r**:
    Growth curve line graph.
    Compares the growth rate of intercept and lag models.
    Overlays the average biomass per age from the field data scatterplot.
    Inputs:


- **3_future_predictions.r**:
    Barplot 1: Compares the biomass gain by 2050 for:
        - random 5% of pastureland
        - 5% with top regrowth potential
        - all secondary forests
    Barplot 2: Shows current area of:
        - 5% of pastureland
        - secondary forests
    Barplot 3: Shows current biomass stocked in:
        - random 5% of pastureland
        - secondary forests
    Shapefile 1: Predicted biomass gain by 2050 for all pastureland.
    Shapefile 2: Predicted biomass gain by 2050 for all secondary forests.
    Inputs:
      - "grid_1k_amazon_secondary": all CSVs in directory
      - "grid_1k_amazon_pastureland": all CSVs in directory


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