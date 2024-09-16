# scripts/nelder_mead_fit.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, cross_validate, GridSearchCV, learning_curve
from sklearn.inspection import permutation_importance
from functools import partial
from xgboost import XGBRegressor

from data_utils import load_and_preprocess_data
from model_utils import objective_function
from tuners import optimize_with_ray_tune, optimize_with_skopt




    # if run_nelder_mead:

    #     nm_mse, nm_r2 = evaluate_nelder_mead(X_train, X_test, y_train, y_test, A_train, A_test, initial_params)
    # fit_params = minimize(lambda pars: objective_function(X_train_scaled, y_train, A_train, pars), initial_params, method='Nelder-Mead')

def create_and_evaluate_model(X, y, model, unseen_data, param_grid=None):

    pipeline = Pipeline([
        ('scaler', MinMaxScaler()),
        ('regressor', model)
    ])

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    splits = list(kf.split(X))
    
    if param_grid:
        param_grid = {'regressor__' + k: v for k, v in param_grid.items()}
        grid_search = GridSearchCV(pipeline, param_grid, cv=kf, 
                            scoring='neg_mean_squared_error', refit = False)
        grid_search.fit(X, y)
        cv_results = grid_search.cv_results_
        mse_scores = -cv_results['mean_test_score']
        best_index = grid_search.best_index_

        # Create and fit the best model using the best parameters
        best_params = grid_search.best_params_
        best_model = Pipeline([
            ('scaler', MinMaxScaler()),
            ('regressor', model.set_params(**{k.split('__')[1]: v for k, v in best_params.items()}))
        ])
    else:
        cv_results = cross_validate(pipeline, X, y, cv=kf, 
                            scoring=['neg_mean_squared_error'],
                            return_estimator=True)
        mse_scores = -cv_results['test_neg_mean_squared_error']
        best_index = np.argmin(-cv_results['test_neg_mean_squared_error'])
        best_model = cv_results['estimator'][best_index]
        train_index, test_index = splits[best_index]
    
    # Perform permutation importance
    _, test_indices = list(kf.split(X))[best_index]
    X_test = X.iloc[test_indices]
    y_test = y[test_indices]
    perm_importance = permutation_importance(best_model, X_test, y_test, n_repeats=10, random_state=42)


    # Calculate MSE of the best model in the unseen dataset
    y_pred = best_model.predict(unseen_data.X)
    unseen_r2 = r2_score(unseen_data.y, y_pred)

    # Calculate R2 from MSE
    r2_scores = 1 - (mse_scores / np.var(y))
    mean_r2 = np.mean(r2_scores)
    std_r2 = np.std(r2_scores)

    return best_model, mean_r2, std_r2, unseen_r2, perm_importance




def print_feature_importance(perm_importance, feature_names):
    feature_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': perm_importance.importances_mean,
        'std': perm_importance.importances_std
    })
    feature_importance = feature_importance.sort_values('importance', ascending=False)
    print("\nFeature Importance:")
    print(feature_importance.to_string(index=False))


def main():
    pars = [
        "nearest_mature", "age", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
        "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
        "cwd"
    ]

    X, y, initial_params, _, unseen_data = load_and_preprocess_data("./0_data/non_aggregated_100k_all.csv", pars)


    models = {
        "Linear Regression": LinearRegression(),
        "XGBoost": XGBRegressor(random_state=42),
        "Random Forest": RandomForestRegressor(random_state=42)
    }

    param_grids = {
        "Linear Regression": None,
        "XGBoost": {
            'learning_rate': [0.001, 0.01, 0.1],
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

    for name, model in models.items():
        _, mean_r2, std_r2, final_r2, perm_importance = create_and_evaluate_model(
            X, y, model, unseen_data, param_grids[name]
        )
        
        print(f"\n{name} Results:")
        print(f"Cross-validation R2: {mean_r2:.3f} (±{std_r2:.3f})")
        print(f"Final model R2: {final_r2:.3f}")
        
        print_feature_importance(perm_importance, X.columns)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # nelder_mead_pars = [par for par in pars if par != "nearest_mature"]

    # X, y, initial_params, A = load_and_preprocess_data("./3_python_scripts/non_aggregated_5000.csv", nelder_mead_pars)

    # # Create a partial function for the objective function
    # partial_objective = partial(objective_function, X=X, y=y, A=A)

    # # Define optimizers
    # optimizers = {
    #     "Initial (Non-tuned)": lambda: initial_params,
    #     "Ray Tune": lambda: optimize_with_ray_tune(partial_objective, nelder_mead_pars),
    #     "Skopt": lambda: optimize_with_skopt(partial_objective, nelder_mead_pars)
    # }

    # # Run optimizers and collect results
    # results = {}
    # for name, optimizer in optimizers.items():
    #     print(f"\nOptimizing with {name}...")
    #     best_params = optimizer()
    #     print(f"Best parameters found by {name}:", best_params)

    #     if name == "Initial (Non-tuned)":
    #         params = initial_params
    #     else:
    #         # Convert best_params to the format your optimization function expects
    #         params = np.array([
    #             best_params["B0"],
    #             best_params["theta"],
    #             *[best_params[f"coeff_{i}"] for i in range(len(nelder_mead_pars))]
    #         ])

    #     # Perform cross-validation
    #     cv_results, _ = cross_validate(X, y, params, A, run_nelder_mead=True)
    #     results[name] = cv_results

    # # Print comparison of results
    # print("\nComparison of Results:")
    # for name, result in results.items():
    #     print(f"\n{name} Optimizer:")
    #     print(f"Cross-validation MSE: {result['nelder_mead_mse_mean']:.4f} (+/- {result['nelder_mead_mse_std']:.4f})")
    #     print(f"Cross-validation R^2: {result['nelder_mead_r2_mean']:.4f} (+/- {result['nelder_mead_r2_std']:.4f})")


if __name__ == "__main__":
    main()