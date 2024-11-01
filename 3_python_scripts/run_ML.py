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
        remove_landscape: bool = False
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
    
    # Convert 'topography' and 'ecoreg' to categorical if present
    for col in ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']:
        if col in df.columns:
            df[col] = df[col].astype('category')

    if pars is None:
        pars = df.columns.tolist()
        if not keep_all_data:
            pars = [col for col in pars if col not in ["biome", "biomass", "latitude", "longitude"]]

    # Remove land use related parameters if the switch is on
    if remove_land_use:
        keywords = ['lulc', 'LU', 'fallow', 'num_fires']
        pars = [par for par in pars if not any(keyword in par for keyword in keywords)]
    # Remove land use related parameters if the switch is on
    if remove_landscape:
        keywords = ['nearest_mature_biomass', 'sur_cover', 'distance']
        pars = [par for par in pars if not any(keyword in par for keyword in keywords)]

    # # Identify columns with fewer than 200 non-zero values
    # low_value_columns = [col for col in df.columns if (df[col] != 0).sum() < 200]

    # # Create a temporary combined category column for adjusted sampling
    # df['stratify_col'] = df[low_value_columns].astype(str).agg('-'.join, axis=1)
    # class_counts = df['stratify_col'].value_counts()
    # df['sampling_weight'] = df['stratify_col'].map(1 / class_counts)
    # df_sample = df.sample(n=10000, weights='sampling_weight', random_state=42, replace=False)
    # df_sample = df_sample.drop(columns=['stratify_col', 'sampling_weight'])

    # df = df.sample(10000, random_state = 42)

    X = df[pars]
    y = df['biomass'].values

    return X, y



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

    # Placeholder for predictions (initialized with NaNs to identify unfilled indices)
    predictions = np.empty(len(X))
    predictions[:] = np.nan

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
    mean_r2 = np.mean(r2_scores)
    std_r2 = np.std(r2_scores)

    # Fit the model on the entire dataset
    final_model = Pipeline([
        ('scaler', MinMaxScaler()),
        ('regressor', model.named_steps['regressor'])
    ])

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    final_model.fit(X_train, y_train)
    final_r2 = r2_score(y_test, final_model.predict(X_test))

    return mean_r2, std_r2, final_r2, perm_importance

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


def run_model(datasource, models, param_grids, biome = 4, remove_land_use = False, remove_landscape = False):
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
    X, y = load_and_preprocess_data(filepath, \
                        biome = biome,
                        final_sample_size = 8000,
                        remove_land_use = remove_land_use,
                        remove_landscape = remove_landscape)

    rows = pd.DataFrame()
    # Perform regression analysis for each model
    for name, model in models.items():
        mean_r2, std_r2, final_r2, perm_importance = regression_cv(
            X, y, model, param_grids[name]
        )
        
        print_feature_importance(perm_importance, X.columns)

        row = {
            'datasource': datasource,
            'biome': biome,
            'model': name,
            'mean_r2': mean_r2,
            'std_r2': std_r2,
            'final_r2': final_r2
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
    # intervals = ["all", "5yr", "10yr", "15yr"]
    # biomes = [4]
    # intervals = ["all"]
    # types = ["aggregated", "non_aggregated"]
    # results_df = pd.DataFrame()
    # for biome in biomes:
    #     for interval in intervals:
    #         for data_type in types:
    #             datasource = f"{data_type}_{interval}"
    #             results = run_model(datasource, 
    #                                 {"XGBoost": models["XGBoost"]},
    #                                 {"XGBoost": param_grids["XGBoost"]},
    #                                  biome)
    #             results_df = pd.concat([results_df, results], ignore_index=True)

    # results_df.to_csv("./0_results/biome_interval_aggregated.csv", index=False)
    # print("Saved to CSV.")

    # datasources = ["aggregated_all", "ESA_fire"]
    # results_df = pd.DataFrame()
    # for datasource in datasources:
    #     results = run_model(datasource, 
    #                                 {"XGBoost": models["XGBoost"]},
    #                                 {"XGBoost": param_grids["XGBoost"]},
    #                                 biome = 1,
    #                                 remove_land_use = True
    #                                 )
    #     results_df = pd.concat([results_df, results], ignore_index=True)

    # # results_df.to_csv("./0_results/mapbiomas_eu.csv", index=False)
    # # print("Saved to CSV.")



    lulc_presence = [True, False]
    landscape_presence = [True, False]
    results_df = pd.DataFrame()
    for remove_landscape in landscape_presence:
        for remove_land_use in lulc_presence:
            results = run_model("aggregated_all", 
                                        {"XGBoost": models["XGBoost"]},
                                        {"XGBoost": param_grids["XGBoost"]},
                                        biome = 1,
                                        remove_land_use = remove_land_use,
                                        remove_landscape = remove_landscape
                                        )
            results_df = pd.concat([results_df, results], ignore_index=True)

    # results_df.to_csv("./0_results/mapbiomas_eu.csv", index=False)
    # print("Saved to CSV.")



    # filepaths = ["mapbiomas_fire", "ESA_fire"]
    # results_df = pd.DataFrame()
    # for datasource in datasources:
    #     results = run_model(datasource, models[2], param_grids[2])
    #     results_df = pd.concat([results_df, results], ignore_index=True)

    # results_df.to_csv("./0_results/mapbiomas_eu.csv", index=False)
    # print("Saved to CSV.")