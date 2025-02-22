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
│   ├── gee_6_write_csv.ipynb
│   └── test_scripts
│       ├── gee_histogram.ipynb
│       ├── gee_manage_assets.py
│       ├── gee_visualize_Planet.ipynb
│       └── modis.py
├── 2_R_scripts
│   ├── 1_data_utils.r
│   ├── 2_model_utils.r
│   └── 3_run.r
├── 3_python_scripts
│   ├── data_utils.py
│   ├── model_utils.py
│   ├── run.py
│   └── tuners.py
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

- **gee_1_secondary.ipynb**:
    - imports secondary forest age data from TMF and Mapbiomas
    - removes pixels with ages that don't match the IPCC estimates
    - removes isolated pixels (keeps only)
    - imports and reprojects GEDI L2A, GEDI L4A and ESA CCI Biomass data

- **gee_2_categorical.ipynb**:
    - imports:
        - Indigenous Lands
        - Protected Areas
        - Ecoregions
        - Biomes
    - exports "categorical" to GEE assets
    - imports:
      - 

- **gee_3_climate.ipynb**:
- **gee_4_mature.ipynb**:
- **gee_5_land_use.ipynb**:
- **gee_6_write_csv.ipynb**:
- **test_scripts/**: 
  - **gee_histogram.ipynb**:
  - **gee_manage_assets.py**:
  - **gee_visualize_Planet.ipynb**:
  - **modis.py**:


### 2_R_scripts/
- **1_data_utils.r**:
- **2_model_utils.r**:
- **3_run.r**:

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

3. Run the data processing script:
    ```sh
    python 3_python_scripts/data_utils.py
    ```

4. Train the models:
    ```sh
    python 3_python_scripts/run.py
    ```

5. Evaluate the models:
    ```sh
    python 3_python_scripts/model_eval.py
    ```

## Dependencies
The project requires the following Python packages:
- pandas
- numpy
- scikit-learn
- matplotlib
- seaborn
- google-earth-engine

Refer to `requirements.txt` for the complete list of dependencies.

## Contributing
Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.
