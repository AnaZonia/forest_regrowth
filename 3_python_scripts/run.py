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
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from functools import partial

from data_utils import load_and_preprocess_data, make_initial_parameters, format_best_params
from model_utils import regression_cv, cross_validate_nelder_mead, nelder_mead_lag, nelder_mead_B0_theta
from tuners import optimize_with_grid_search, optimize_with_ray_tune, optimize_with_skopt


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


def regression_main():
    """
    Main function to perform regression analysis using different models.
    Loads data, performs 5-fold cross-validation, and prints results.
    """
    pars = ["age", "cwd"]

    X, y, _, _, unseen_data = load_and_preprocess_data("./0_data/non_aggregated_100k_all.csv", pars)

    models = {
        "Linear Regression": LinearRegression(),
        "XGBoost": XGBRegressor(random_state = 42),
        "Random Forest": RandomForestRegressor(random_state = 42),
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
    
    figures = []
    for name, model in models.items():
        mean_r2, std_r2, unseen_r2, perm_importance, fig = regression_cv(
            X, y, model, unseen_data, name, param_grids[name]
        )
        
        print(f"\n{name} Results:")
        print(f"Cross-validation R2: {mean_r2:.3f} (±{std_r2:.3f})")
        print(f"Unseen data R2: {unseen_r2:.3f}")
        
        print_feature_importance(perm_importance, X.columns)
        figures.append((name, fig))

    # Display all figures at the end
    for name, fig in figures:
        fig.suptitle(name)
        plt.show()



def nelder_mead_main(tune = False, func_form = 'lag'):
    """
    Main function to perform Nelder-Mead optimization with optional parameter tuning.

    Parameters:
    - tune: bool, whether to perform parameter tuning using different methods.
    - func_form: str, the functional form to use ('lag' or 'B0_theta').

    Returns:
    - None, but prints the results of cross-validation and unseen data R2 scores.
    """
    pars = ["age", "cwd"]

    X, y, A, unseen_data = load_and_preprocess_data("./0_data/non_aggregated_100k_all.csv", pars)
    initial_params = make_initial_parameters(pars, y, lag = True)

    # Define tuners
    init_param_tuners = {
        "Initial (Non-tuned)": lambda: initial_params
    }

    if tune:
        # Create a partial function for the objective function
        model = partial(nelder_mead_lag, X = X, y = y)
        # Add tuning methods
        init_param_tuners.update({
            "Grid Search": lambda: optimize_with_grid_search(model, initial_params),
            "Ray Tune": lambda: optimize_with_ray_tune(model, initial_params),
            "Skopt": lambda: optimize_with_skopt(model, initial_params)
        })


    # Run tuners and collect results
    figures = []
    for name, init_params in init_param_tuners.items():
        print(f"\nTuning with {name}...")
        best_params = init_params()
        print(f"Best parameters found by {name}:", best_params)

        if name == "Initial (Non-tuned)":
            params = initial_params
        else:
            # Convert best_params to the format your optimization function expects
            params = format_best_params(best_params, pars, func_form)

        if func_form == "B0_theta":
            nelder_mead_func = nelder_mead_B0_theta
        elif func_form == "lag":
            nelder_mead_func = nelder_mead_lag

        mean_score, std_score, unseen_r2, fig = cross_validate_nelder_mead(
            X, y, A, params, unseen_data, name, nelder_mead_func
        )
        
        print(f"\n{name} Results:")
        print(f"Cross-validation values: {mean_score:.3f} (±{std_score:.3f})")
        print(f"Unseen data R2: {unseen_r2:.3f}")
        
        figures.append((name, fig))

if __name__ == "__main__":
    nelder_mead_main()