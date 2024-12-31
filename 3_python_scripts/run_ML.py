"""
Nelder-Mead Optimization and Regression Analysis Script

This script performs regression analysis and Nelder-Mead optimization on a dataset.
It includes functionality for:
1. Standard regression analysis using Linear Regression, XGBoost, and Random Forest
2. Nelder-Mead optimization with optional hyperparameter tuning
3. Visualization of results

Modules:
- data_utils: Contains utility functions for loading, preprocessing, and handling data.
- model_utils: Contains functions for performing regression and Nelder-Mead optimization.
- tuners: Contains functions for hyperparameter tuning using various methods.

Functions:
- print_feature_importance: Prints the feature importance based on permutation importance.
- regression_main: Main function to perform regression analysis using different models.
- nelder_mead_main: Main function to perform Nelder-Mead optimization with optional parameter tuning.

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, cross_validate, GridSearchCV, learning_curve
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit

def load_and_preprocess_data(
        filepath: str, 
        pars = None, 
        biome = "both",
        keep_all_data: bool = False,
        final_sample_size: int = 10000,
        remove_land_use: bool = False,
        remove_landscape: bool = False,
        keep_only_land_use_and_landscape: bool = False
    ) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Load and preprocess data from a CSV file.

    Args:
        filepath (str): Path to the CSV file.
        pars (List[str]): List of parameter names to keep.
        keep_all_data (bool): Flag to keep all data or subset by 'biome'. Defaults to False.
        final_sample_size (int): Total sample size to use after stratified sampling. Defaults to 10000.

    Returns:
        Tuple[pd.DataFrame, np.ndarray]: Features (X) and target (y).
    """
    df = pd.read_csv(filepath)        
    
    if biome != "both":
        df = df[df['biome'] == biome]
    # if biome == 4:
    #     df = df.drop(columns = ['mean_aet', 'cwd']) # multicollinearity
    # if biome == "both":
    #     df = df.drop(columns = ['mean_aet']) # multicollinearity
    df = df.loc[:, ~df.columns.str.contains(r"\d{4}")]
    
    # Identify columns with fewer than 200 non-zero values
    low_value_columns = [col for col in df.columns if (df[col] != 0).sum() < 100]
    df = df.drop(columns = low_value_columns)
    # df = df.sample(10000, random_state = 42)
    # # Create a temporary combined category column for adjusted sampling
    # df['stratify_col'] = df[low_value_columns].astype(str).agg('-'.join, axis=1)
    # class_counts = df['stratify_col'].value_counts()
    # df['sampling_weight'] = df['stratify_col'].map(1 / class_counts)
    # df_sample = df.sample(n=10000, weights='sampling_weight', random_state=42, replace=False)
    # df_sample = df_sample.drop(columns=['stratify_col', 'sampling_weight'])

    # Convert 'topography' and 'ecoreg' to categorical if present
    for col in ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']:
        if col in df.columns:
            df[col] = df[col].astype('category')
            # Filter out categories with fewer than 50 rows
            counts = df[col].value_counts()
            valid_categories = counts[counts >= 50].index
            df = df[df[col].isin(valid_categories)]

    if pars is None:
        pars = df.columns.tolist()
        if not keep_all_data:
            pars = [col for col in pars if col not in ["biome", "biomass", "latitude", "longitude"]]

    # Define keywords for land use and landscape parameters, and include "age"
    land_use_keywords = ['lulc', 'LU', 'fallow', 'num_fires']
    landscape_keywords = ['sur_cover', 'distance', 'nearest_mature_biomass']

    essential_keywords = land_use_keywords + landscape_keywords + ['age']

    # Filter parameters based on the switches
    if remove_land_use:
        pars = [par for par in pars if not any(keyword in par for keyword in land_use_keywords)]
        
    if remove_landscape:
        pars = [par for par in pars if not any(keyword in par for keyword in landscape_keywords)]

    # New switch to keep only land use, landscape, and "age" parameters
    if keep_only_land_use_and_landscape:
        pars = [par for par in pars if any(keyword in par for keyword in essential_keywords)]

    X = df[pars]
    y = df['biomass'].values
    lat_lon = df[['latitude', 'longitude']]

    return X, y, lat_lon



def plot_learning_curves(X, y, model, name, cv = 5): 
    """
    Plot learning curves for a given model pipeline.
    
    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        model: The model to train.
        name (str): Name of the model for the plot title.
        cv (int): Number of cross-validation folds.
    """
    train_sizes, train_scores, test_scores = learning_curve(
        model, X, y, cv = cv, scoring = 'neg_mean_squared_error',
        train_sizes = np.linspace(0.1, 1.0, 10), n_jobs = -1)
    
    # Calculate mean and standard deviation for training set scores
    train_mean = -np.mean(train_scores, axis = 1)
    train_std = np.std(train_scores, axis = 1)

    # Calculate mean and standard deviation for test set scores
    test_mean = -np.mean(test_scores, axis = 1)
    test_std = np.std(test_scores, axis = 1)

    # Plot learning curve
    plt.figure(figsize = (10, 6))
    plt.fill_between(train_sizes, train_mean - train_std, train_mean + train_std, alpha = 0.1, color = "r")
    plt.fill_between(train_sizes, test_mean - test_std, test_mean + test_std, alpha = 0.1, color = "g")
    plt.plot(train_sizes, train_mean, 'o-', color = "r", label = "Training score")
    plt.plot(train_sizes, test_mean, 'o-', color = "g", label = "Cross-validation score")
    
    plt.xlabel("Training examples")
    plt.ylabel("Mean Squared Error")
    plt.title(f"Learning Curves for {name}")
    plt.legend(loc = "best")
    plt.grid(True)
    
    return plt



def regression_cv(X, y, model, param_grid = None):
    """
    Perform cross-validation for regression models.

    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        model: The regression model to train.
        param_grid (dict): Parameter grid for hyperparameter tuning.

    Returns:
        tuple: Mean R2, standard deviation of R2, permutation importance, predictions array.
    """
    model = Pipeline([
        ('scaler', MinMaxScaler()),
        ('regressor', model)
    ])

    kf = KFold(n_splits = 5, shuffle = True, random_state = 42)
    splits = list(kf.split(X))
    
    if param_grid:
        param_grid = {'regressor__' + k: v for k, v in param_grid.items()}
        grid_search = GridSearchCV(model, param_grid, cv = kf, 
                            scoring = 'neg_mean_squared_error', refit = True)
        grid_search.fit(X, y)

        best_params = grid_search.best_params_
        model = model.set_params(**best_params)
    else:
        model = model.fit(X, y)


    # After performing grid search, perform cross-validation
    cv_results = cross_validate(model, X, y, cv = kf, 
                        scoring = ['neg_mean_squared_error'],
                        return_estimator = True)
    mse_scores = -cv_results['test_neg_mean_squared_error']

    # Perform permutation importance
    best_index = np.argmin(-cv_results['test_neg_mean_squared_error'])
    _, test_indices = splits[best_index]
    X_test = X.iloc[test_indices]
    y_test = y[test_indices]
    best_model = cv_results['estimator'][best_index]
    perm_importance = permutation_importance(best_model, X_test, y_test, n_repeats = 10, random_state = 42)

    # Calculate R2 from MSE
    r2_scores = 1 - (mse_scores / np.var(y))
    r2_mean = np.mean(r2_scores)
    r2_sd = np.std(r2_scores)

    # Gather predictions for each fold in the original order of X
    # Placeholder for predictions (initialized with NaNs to identify unfilled indices)

    # Fit the model on the entire dataset
    # final_model = Pipeline([
    #     ('scaler', MinMaxScaler()),
    #     ('regressor', model.named_steps['regressor'])
    # ])
    # final_model.fit(X, y)

    predictions = np.empty(len(X))
    predictions[:] = np.nan
    for i, estimator in enumerate(cv_results['estimator']):
        _, test_indices = list(kf.split(X))[i]
        fold_predictions = estimator.predict(X.iloc[test_indices])
        predictions[test_indices] = fold_predictions

    # Calculate R-squared between all biomass values (y) and all predictions across folds
    r2_final = r2_score(y, predictions)

    return r2_mean, r2_sd, r2_final, perm_importance, predictions

def print_feature_importance(perm_importance, feature_names):
    """
    Print the feature importance based on permutation importance.

    Parameters:
    - perm_importance: PermutationImportance object containing importances.
    - feature_names: List of feature names.
    """
    feature_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': perm_importance.importances_mean,
        'std': perm_importance.importances_std
    })
    feature_importance = feature_importance.sort_values('importance', ascending = False)
    print("\nFeature Importance:")
    print(feature_importance.to_string(index = False))


def run_model(datasource, models, param_grids, biome = 4, remove_land_use = False, remove_landscape = False, keep_only_land_use_and_landscape = False):
    """
    Run regression models on a specified datasource with configurable biome and land-use parameter removal.

    Args:
        datasource (str): The name of the data source file (without ".csv").
        biome (int): The biome to filter the data on. Defaults to 1.
        remove_land_use (bool): Whether to remove land-use related parameters. Defaults to False.

    Returns:
        pd.DataFrame: DataFrame containing the results for each model.
    """

    filepath = f"./0_data/{datasource}.csv"

    # Load and preprocess data for the given biome
    X, y, _ = load_and_preprocess_data(filepath, \
                        biome = biome,
                        final_sample_size = 8000,
                        remove_land_use = remove_land_use,
                        remove_landscape = remove_landscape,
                        keep_only_land_use_and_landscape = keep_only_land_use_and_landscape)

    rows = pd.DataFrame()
    # Perform regression analysis for each model
    for name, model in models.items():
        r2_mean, r2_sd, r2_final, perm_importance, _ = regression_cv(
            X, y, model, param_grids[name]
        )
        
        print_feature_importance(perm_importance, X.columns)

        row = {
            'datasource': datasource,
            'biome': biome,
            'model': name,
            'land_use': not remove_land_use,
            'landscape': not remove_landscape,
            'keep_only_land_use_and_landscape': keep_only_land_use_and_landscape,
            'r2_mean': r2_mean,
            'r2_sd': r2_sd,
            'r2_final': r2_final
        }

        print(row)
        # Append result to the DataFrame
        rows = pd.concat([rows, pd.DataFrame([row])], ignore_index=True)

    return rows



if __name__ == "__main__":

    models = {
        "Linear Regression": LinearRegression(),
        "XGBoost": XGBRegressor(random_state = 42),
        "Random Forest": RandomForestRegressor(random_state = 42)
    }

    param_grids = {
        "Linear Regression": None,
        "XGBoost": {
            'learning_rate': [0.01, 0.1],
            'n_estimators': [100, 500],
            'max_depth': [3, 7],
            'min_child_weight': [1, 5]
        },
        "Random Forest": {
            'n_estimators': [100, 300],
            'max_depth': [None, 20],
            'min_samples_split': [2, 10],
            'min_samples_leaf': [1, 4]
        }
    }

    # Define biomes, intervals, and types
    # biomes = [1, 4, "both"]
    # types = ["aggregated_all", "aggregated_5yr", "aggregated_10yr", "aggregated_15yr"]
    # types = ["aggregated_all_tilewise", "aggregated_5yr_tilewise", "aggregated_10yr_tilewise", "aggregated_15yr_tilewise"]
    # remove_lulc = [True, False]
    # remove_landscape = [True, False]
    # # Calculate total number of iterations
    # total_iterations = len(biomes) * len(types) * len(landscape_presence) * len(lulc_presence)
    # current_iteration = 0  # Initialize the counter

    # results_df = pd.DataFrame()
    # for biome in biomes:
    #     for data_type in types:
    #         datasource = f"{data_type}"
    #         results = run_model(datasource, 
    #                             {"XGBoost": models["XGBoost"]},
    #                             {"XGBoost": param_grids["XGBoost"]},
    #                             biome,
    #                             keep_only_land_use_and_landscape = True
    #                             )
    #         # Append results and update progress
    #         results_df = pd.concat([results_df, results], ignore_index=True)
    #         # current_iteration += 1
    #         # print(f"Progress: {current_iteration}/{total_iterations} ({(current_iteration / total_iterations) * 100:.2f}%) - Saved to CSV for  biome '{biome}'.")


    #         # Save CSV after each `data_type` and `interval` iteration set
    #         results_df.to_csv("./0_results/keep_only_land_use_and_landscape.csv", index=False)
    #         print("Data saved to CSV.")



    # # Load and preprocess data for the given biome
    # X, y, lat_lon = load_and_preprocess_data("./0_data/non_aggregated_all.csv", 
    #                                         biome=1,
    #                                         final_sample_size=8000,
    #                                         remove_land_use=False,
    #                                         remove_landscape=False)

    # # Run regression with cross-validation
    # _, _, _, _, predictions = regression_cv(
    #     X, y, models["XGBoost"],
    #         param_grids["XGBoost"]
    # )

    # # Convert y and lat_lon to DataFrames for consistency and merging
    # y_df = pd.DataFrame(y, columns=["biomass"])
    # lat_lon_df = pd.DataFrame(lat_lon)

    # # Concatenate all data with predictions into a single DataFrame
    # results_df = pd.concat([pd.DataFrame(X).reset_index(drop=True), 
    #                         y_df.reset_index(drop=True),
    #                         lat_lon_df.reset_index(drop=True),
    #                         pd.DataFrame(predictions, columns=["predictions"]).reset_index(drop=True)], axis=1)

    # # Save the resulting DataFrame as CSV
    # results_df.to_csv("non_aggregated_all_pred.csv", index=False)



    # # Define biomes, intervals, and types
    # biomes = [1]
    # types = ["aggregated"]#, "eu"]
    # intervals = ["all"]
    # remove_lulc = [True, False]
    # remove_landscape = [True, False]
    # keep_only_landscape_lulc = [True, False]
    # # Calculate total number of iterations
    # total_iterations = len(biomes) * len(types) * len(remove_landscape) * len(remove_lulc) * (1 + len(intervals))
    # current_iteration = 0  # Initialize the counter

    # results_df = pd.DataFrame()
    # for biome in biomes:
    #     for data_type in types:
    #         for interval in intervals:
    #             for landscape in remove_landscape:
    #                 for land_use in remove_lulc:
    #                     if not landscape and not land_use:
    #                         for keep_only in keep_only_landscape_lulc:
    #                             datasource = f"{data_type}"
    #                             datasource = f"{data_type}_{interval}"
    #                             results = run_model(datasource, 
    #                                                 {"XGBoost": models["XGBoost"]},
    #                                                 {"XGBoost": param_grids["XGBoost"]},
    #                                                 biome,
    #                                                 remove_land_use = land_use,
    #                                                 remove_landscape = landscape,
    #                                                 keep_only_land_use_and_landscape = keep_only
    #                                                 )
    #                             # Append results and update progress
    #                             results_df = pd.concat([results_df, results], ignore_index=True)
    #                             current_iteration += 1
    #                             print(f"Progress: {current_iteration}/{total_iterations} ({(current_iteration / total_iterations) * 100:.2f}%) - Saved to CSV for data type '{data_type}', interval '{interval}', biome '{biome}'.")

    #                             # Save CSV after each `data_type` and `interval` iteration set
    #                             results_df.to_csv("./0_results/all_iters.csv", index=False)
    #                             print("Data saved to CSV.")
    #                     else:
    #                         datasource = f"{data_type}"
    #                         datasource = f"{data_type}_{interval}"
    #                         results = run_model(datasource, 
    #                                             {"XGBoost": models["XGBoost"]},
    #                                             {"XGBoost": param_grids["XGBoost"]},
    #                                             biome,
    #                                             remove_land_use = land_use,
    #                                             remove_landscape = landscape
    #                                             )
    #                         # Append results and update progress
    #                         results_df = pd.concat([results_df, results], ignore_index=True)
    #                         current_iteration += 1
    #                         print(f"Progress: {current_iteration}/{total_iterations} ({(current_iteration / total_iterations) * 100:.2f}%) - Saved to CSV for data type '{data_type}', interval '{interval}', biome '{biome}'.")

    #                         # Save CSV after each `data_type` and `interval` iteration set
    #                         results_df.to_csv("./0_results/all_iters_amaz.csv", index=False)
    #                         print("Data saved to CSV.")




    # # Define biomes, intervals, and types
    # biomes = [1]
    # types = ["aggregated_all"]

    # results_df = pd.DataFrame()
    # for biome in biomes:
    #     for data_type in types:
    #         datasource = f"{data_type}"
    #         results = run_model(datasource, 
    #                             {"XGBoost": models["XGBoost"]},
    #                             {"XGBoost": param_grids["XGBoost"]},
    #                             biome
    #                             # remove_land_use = land_use,
    #                             # remove_landscape = landscape
    #                             )
    #         # Append results and update progress
    #         results_df = pd.concat([results_df, results], ignore_index=True)
    #         # Save CSV after each `data_type` and `interval` iteration set
    #         results_df.to_csv("./0_results/all_iters.csv", index=False)
    #         print("Data saved to CSV.")



# Define biomes, intervals, and types
biomes = [1]
types = ["aggregated"]  # ["aggregated", "eu"] can be expanded if needed
intervals = ["all"]
# Only the desired configurations
configurations = [
    {"keep_only_land_use_and_landscape": True, "remove_landscape": False, "remove_land_use": False},
    {"keep_only_land_use_and_landscape": False, "remove_landscape": True, "remove_land_use": True},
]

# Calculate total number of iterations
total_iterations = len(biomes) * len(types) * len(intervals) * len(configurations)
current_iteration = 0  # Initialize the counter

results_df = pd.DataFrame()
for biome in biomes:
    for data_type in types:
        for interval in intervals:
            for config in configurations:
                # Set configuration parameters
                keep_only = config["keep_only_land_use_and_landscape"]
                landscape = config["remove_landscape"]
                land_use = config["remove_land_use"]

                datasource = f"{data_type}_{interval}"
                results = run_model(
                    datasource,
                    {"XGBoost": models["XGBoost"]},
                    {"XGBoost": param_grids["XGBoost"]},
                    biome,
                    remove_land_use=land_use,
                    remove_landscape=landscape,
                    keep_only_land_use_and_landscape=keep_only
                )

                # Append results and update progress
                results_df = pd.concat([results_df, results], ignore_index=True)
                current_iteration += 1
                print(f"Progress: {current_iteration}/{total_iterations} "
                      f"({(current_iteration / total_iterations) * 100:.2f}%) - "
                      f"Saved to CSV for data type '{data_type}', interval '{interval}', biome '{biome}'.")

                # Save CSV after each iteration
                results_df.to_csv("./0_results/all_iters_2.csv", index=False)
                print("Data saved to CSV.")
