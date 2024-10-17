## Project Overview
CapoERA (Computational Algorithm for Post-agricultural Ecosystem Regrowth Analysis) is a model that predicts the age of secondary forests in Brazil based on satellite data. This script imports and processes remote sensing data from Google Earth Engine, and compares different models to predict the biomass of secondary forests based on climate and previous land use.


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
├── 4_stan_scripts
│   ├── model.stan
│   └── run_model_stan.R
├── 5_testing_scripts
│   ├── conceptualize_growth_curves.py
│   ├── eu_mapbiomas_exploration.py
│   ├── hyper_tune_bayes_skopt.py
│   ├── investigate_mature_biomass.qmd
│   ├── multicollinearity.py
│   └── plot_histogram.py
├── methods.txt
├── README.md
└── requirements.txt


### 1_gee_scripts/
- **gee_0_utils.py**: Utility functions for Google Earth Engine (GEE) scripts.
- **gee_1_secondary.ipynb**: Notebook for processing secondary forest data.
- **gee_2_categorical.ipynb**: Notebook for processing categorical data.
- **gee_3_climate.ipynb**: Notebook for processing climate data.
- **gee_4_mature.ipynb**: Notebook for processing mature forest data.
- **gee_5_land_use.ipynb**: Notebook for processing land use data.
- **gee_6_write_csv.ipynb**: Notebook for writing data to CSV files.
- **test_scripts/**: Contains test scripts for various GEE functionalities.
  - **gee_histogram.ipynb**: Notebook for creating histograms.
  - **gee_manage_assets.py**: Script for managing GEE assets.
  - **gee_visualize_Planet.ipynb**: Notebook for visualizing Planet data.
  - **modis.py**: Script for processing MODIS data.

### 2_R_scripts/
- **1_data_utils.r**: Utility functions for data processing in R.
- **2_model_utils.r**: Utility functions for model training in R.
- **3_run.r**: Script to run the R models.

### 3_python_scripts/
- **data_utils.py**: Utility functions for data processing in Python.
- **model_utils.py**: Utility functions for model training in Python.
- **run.py**: Script to run the Python models.
- **tuners.py**: Script for hyperparameter tuning.

### 4_stan_scripts/
- **model.stan**: Stan model definition.
- **run_model_stan.R**: Script to run the Stan model using R.

### 5_testing_scripts/
- **conceptualize_growth_curves.py**: Script to conceptualize growth curves.
- **eu_mapbiomas_exploration.py**: Script for exploring EU MapBiomas data.
- **hyper_tune_bayes_skopt.py**: Script for hyperparameter tuning using Bayesian optimization.
- **investigate_mature_biomass.qmd**: Quarto document for investigating mature biomass.
- **multicollinearity.py**: Script to check for multicollinearity.
- **plot_histogram.py**: Script to plot histograms.

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
