## Project Overview
CapoERA (Computational Algorithm for Post-agricultural Ecosystem Regrowth Analysis) is a model that predicts the age of secondary forests in Brazil based on satellite data. This script imports and processes remote sensing data from Google Earth Engine, and compares different models to predict the biomass of secondary forests based on climate and previous land use.

```
forest_regrowth
├── 1_gee
│   ├── gee_0_utils.py
│   ├── gee_1_secondary.ipynb
│   ├── gee_2_categorical.ipynb
│   ├── gee_3_climate.ipynb
│   ├── gee_4_mature.ipynb
│   ├── gee_5_land_use.ipynb
│   └── gee_6_write_csv.ipynb
├── 2_R_scripts
│   ├── fit_1_import_data.r
│   ├── fit_2_functions.r
│   └── fit_3_run_model.r
├── 3_python_scripts
│   ├── data_utils.py
│   ├── model_utils.py
│   ├── run.py
│   └── tuners.py
├── 4_visualization
│   ├── 1_visualize_data.ipynb
|
├── methods.txt
├── README.md
└── requirements.txt
```

### 1_gee_scripts/
- **gee_0_utils.py**:
    - defines project date range (1985-2020) and imports region of interest
    - defines export_image function
    - imports files for data processing (from Google Earth Engine)
    - select only pixels with exclusively the desired land use histories (exclude all instances of land use types we are not interested in)

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
    - reprojects and exports GEDI L2A, GEDI L4A and ESA CCI Biomass data

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
    - IPCC future climate trajectories

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

### 2_R_scripts/
- **0_multicollinearity.r**:
- **fit_1_import_data.r**:
  - Imports 
- **fit_2_functions.r**:
- **fit_3_run_model.r**:

### 3_python_scripts/
- **data_utils.py**:
- **model_utils.py**:
- **run.py**:
- **tuners.py**:

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
    ```R
    renv::restore()
    ```

Refer to `requirements.txt` for the complete list of dependencies.


