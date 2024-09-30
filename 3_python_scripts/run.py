# scripts/nelder_mead_fit.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from functools import partial

from data_utils import load_and_preprocess_data
from model_utils import regression_cv, nelder_mead_cv, nelder_mead
from tuners import optimize_with_ray_tune, optimize_with_skopt, optimize_with_grid_search


def print_feature_importance(perm_importance, feature_names):
    feature_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': perm_importance.importances_mean,
        'std': perm_importance.importances_std
    })
    feature_importance = feature_importance.sort_values('importance', ascending=False)
    print("\nFeature Importance:")
    print(feature_importance.to_string(index=False))


def regression_main():
    pars = [
        "age", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
        "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
        "cwd"
    ]

    X, y, _, _, unseen_data = load_and_preprocess_data("./0_data/non_aggregated_100k_all.csv", pars)

    models = {
        "Linear Regression": LinearRegression(),
        "XGBoost": XGBRegressor(random_state=42),
        "Random Forest": RandomForestRegressor(random_state=42),

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



def nelder_mead_main():
    pars = [
        "age", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
        "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
        "cwd"
    ]

    X, y, initial_params, A, unseen_data = load_and_preprocess_data("./0_data/non_aggregated_100k_all.csv", pars)

    # Create a partial function for the objective function
    model = partial(nelder_mead, X=X, y=y, A=A)

    # Define optimizers
    optimizers = {
        # "Initial (Non-tuned)": lambda: initial_params,
        "Grid Search": lambda: optimize_with_grid_search(model, initial_params[2:])
        # "Ray Tune": lambda: optimize_with_ray_tune(model, initial_params[2:]),
        # "Skopt": lambda: optimize_with_skopt(model, initial_params[2:])
    }

    # Run optimizers and collect results
    figures = []
    for name, optimizer in optimizers.items():
        print(f"\nOptimizing with {name}...")
        best_params = optimizer()
        print(f"Best parameters found by {name}:", best_params)

        if name == "Initial (Non-tuned)":
            params = initial_params
        else:
            # Convert best_params to the format your optimization function expects
            params = np.array([
                best_params["B0"],
                best_params["theta"],
                *[best_params[f"coeff_{i}"] for i in range(len(pars))]
            ])

        mean_r2, std_r2, unseen_r2 = nelder_mead_cv(
            X, y, A, params, unseen_data, name
        )
        
        print(f"\n{name} Results:")
        print(f"Cross-validation R2: {mean_r2:.3f} (±{std_r2:.3f})")
        print(f"Unseen data R2: {unseen_r2:.3f}")
        
        # figures.append((name, fig))

if __name__ == "__main__":
    regression_main()