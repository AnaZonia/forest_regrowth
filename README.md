## Project Overview
CapoERA (Computational Algorithm for Post-agricultural Ecosystem Regrowth Analysis) is a model that predicts the age of secondary forests in Brazil based on satellite data. This script imports and processes remote sensing data from Google Earth Engine, and compares different models to predict the biomass of secondary forests based on climate and previous land use.

your_project/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── src/
│   ├── __init__.py
│   ├── data_processing.py
│   ├── nelder_mead_model.py
│   ├── random_forest_model.py
│   ├── hyperparameter_tuning.py
│   └── utils.py
│
├── notebooks/
│   ├── exploratory_data_analysis.ipynb
│   └── model_comparison.ipynb
│
├── tests/
│   ├── __init__.py
│   ├── test_data_processing.py
│   ├── test_nelder_mead_model.py
│   └── test_random_forest_model.py
│
├── main.py
├── requirements.txt
└── README.md



### Google Earth Engine Data Processing
- GEE_1_satellite.ipynb - Processes satellite data from MAPBIOMAS and ESA CCI Biomass
	- Fire
	- Land use/land cover
	- Biomass
	- Secondary Forest ages
- GEE_2_categorical.ipynb - Processes categorical data
	- Indigenous Areas
	- Conservation Units
	- Biome/Ecoregion
	- Soil Type
- GEE_3_climate.ipynb - Processes climate data from TerraClim
	- Temperature
	- Precipitation
	- Seasonality
- GEE_4_unifying_data.ipynb - Combines all data into a single dataframe, exported to Google Drive

## Modelling
- fit_1_import_data.r - Imports unified CSV file and prepares climatic variables for modelling
- fit_2_lm_rf_gam.r - Preliminary analyses with linear models, random forests, and generalized additive models	
- fit_3_optim.r - Process-based modelling

## Dependencies
- Google Earth Engine
- Python API
- geemap

install.packages(c("curl", "httr", "openssl", "tidyverse", # required for tidyverse in Fedora
                   "terra", "ggplot2"))
