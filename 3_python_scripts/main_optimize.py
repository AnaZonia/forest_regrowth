# scripts/nelder_mead_fit.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from functools import partial
from xgboost import XGBRegressor
from typing import List, Tuple, Dict
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle

from tuners import optimize_with_ray_tune, optimize_with_skopt


def load_and_preprocess_data(filepath: str, pars: List[str], keep_biome: bool = False, sample_size: int = 10000) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Load and preprocess data from a CSV file.

    Args:
        filepath (str): Path to the CSV file.
        pars (List[str]): List of parameter names to keep.
        sample_size (int, optional): Number of samples to use. Defaults to 10000.

    Returns:
        Tuple[pd.DataFrame, np.ndarray, np.ndarray]: Features (X), target (y), and initial parameters.
    """
    df = pd.read_csv(filepath)

    if not keep_biome:
        df = df[df['biome'] == 4]

    # Step 1: Filter rows where at least one of lulc_sum_41, lulc_sum_39, or num_fires_before_regrowth is non-zero
    condition_df = df[(df['lulc_sum_41'] != 0) | (df['lulc_sum_39'] != 0) | (df['num_fires_before_regrowth'] != 0)]

    # Step 2: Sample at least 2000 rows from this condition_df
    condition_sample = condition_df.sample(n=2000, random_state=42)

    # Step 3: Sample the remaining rows from the rest of the DataFrame
    remaining_rows = sample_size - 2000
    rest_df = df.drop(condition_sample.index)  # Exclude the already sampled rows
    rest_sample = rest_df.sample(n=remaining_rows, random_state=42)

    # Step 4: Combine both samples
    df = pd.concat([condition_sample, rest_sample])

    if keep_biome:
        X = df[pars + ['biome', 'nearest_mature']]
    else:
        X = df[pars]

    y = df['agbd'].values
    A = df['nearest_mature'].values #asymptote
    print((X != 0).sum())

    initial_params = np.zeros(len(pars) + 2)
    initial_params[0] = y.mean()
    initial_params[1] = 1

    return X, y, initial_params, A

def objective_function(X: pd.DataFrame, y: np.ndarray, A: np.ndarray, pars: np.ndarray, mse: bool = True) -> float:
    """
    Calculate the objective function (growth curve + MSE or R^2).

    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        pars (np.ndarray): Model parameters.
        mse (bool, optional): If True, return MSE; otherwise, return R^2. Defaults to True.

    Returns:
        float: MSE or R^2 score.
    """
    B0, theta = pars[:2]
    coeffs = pars[2:]

    adjustment_value = -np.log(1 - (y.mean() / A.mean()))

    k = np.dot(X.values, coeffs) + adjustment_value
    k = np.where(k < 0, adjustment_value, k)

    y_pred = B0 + (A - B0) * (1 - np.exp(-k))**theta

    return mean_squared_error(y, y_pred) if mse else r2_score(y, y_pred)

def optimize_parameters(X: pd.DataFrame, y: np.ndarray, A: np.ndarray, initial_params: np.ndarray) -> np.ndarray:
    """
    Optimize parameters using Nelder-Mead method.

    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        initial_params (np.ndarray): Initial parameter values.

    Returns:
        np.ndarray: Optimized parameters.
    """
    result = minimize(lambda pars: objective_function(X, y, A, pars), initial_params, method='Nelder-Mead')
    return result.x

def scale_features(X_train: pd.DataFrame, X_test: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Scale features using MinMaxScaler, optionally excluding a specified column.

    Args:
        X_train (pd.DataFrame): Training feature matrix.
        X_test (pd.DataFrame): Testing feature matrix.
        exclude_col (str, optional): Column to exclude from scaling. Defaults to None.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Scaled training and testing feature matrices.
    """
    scaler = MinMaxScaler()  

    X_train_for_scaling = X_train
    X_test_for_scaling = X_test

    X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train_for_scaling), columns=X_train_for_scaling.columns, index=X_train.index)
    X_test_scaled = pd.DataFrame(scaler.transform(X_test_for_scaling), columns=X_test_for_scaling.columns, index=X_test.index)
    
    return X_train_scaled, X_test_scaled

def train_and_evaluate_regression(model, X_train: pd.DataFrame, X_test: pd.DataFrame, y_train: np.ndarray, y_test: np.ndarray) -> Tuple[float, float]:
    """
    Train and evaluate a model.

    Args:
        model: The model to train and evaluate.
        X_train (pd.DataFrame): Training feature matrix.
        X_test (pd.DataFrame): Testing feature matrix.
        y_train (np.ndarray): Training target values.
        y_test (np.ndarray): Testing target values.

    Returns:
        Tuple[float, float, object]: MSE, R^2 scores, and trained model.
    """

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    return mse, r2, model

def train_and_evaluate_nelder_mead(X_train: pd.DataFrame, X_test: pd.DataFrame, y_train, y_test, A_train, A_test, initial_params) -> Tuple[float, float]:
    """
    Train and evaluate the Nelder-Mead model.

    Args:
        X_train (pd.DataFrame): Training feature matrix.
        X_test (pd.DataFrame): Testing feature matrix.
        y_train (np.ndarray): Training target values.
        y_test (np.ndarray): Testing target values.
        initial_params (np.ndarray): Initial parameter values.

    Returns:
        Tuple[float, float]: MSE and R^2 scores.
    """
    X_train_scaled, X_test_scaled = scale_features(X_train, X_test)
    fit_params = optimize_parameters(X_train_scaled, y_train, A_train, initial_params)
    mse = objective_function(X_test_scaled, y_test, A_test, fit_params)
    r2 = objective_function(X_test_scaled, y_test, A_test, fit_params, False)
    return mse, r2


def cross_validate(X: pd.DataFrame, y: np.ndarray, initial_params: np.ndarray, A:np.ndarray = None, run_nelder_mead: bool = True) -> Tuple[Dict[str, Dict[str, float]], Dict[str, object]]:
    """
    Perform cross-validation for Nelder-Mead, Linear, XGBoost, and Random Forest models.

    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        initial_params (np.ndarray): Initial parameter values for Nelder-Mead.
        run_nelder_mead (bool): Flag to control whether to run Nelder-Mead optimization.

    Returns:
        Tuple[Dict[str, Dict[str, float]], Dict[str, object]]: Mean and standard deviation of MSE and R^2 for each model, and trained models.
    """
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    
    scores = {
        'nelder_mead': {'mse': [], 'r2': []},
        'linear': {'mse': [], 'r2': []},
        'xgboost': {'mse': [], 'r2': []},
        'random_forest': {'mse': [], 'r2': []}  # Add this line
    }
    
    trained_models = {}

    for train_index, test_index in kf.split(X):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y[train_index], y[test_index]

        if run_nelder_mead:
            # Nelder-Mead model
            A_train, A_test = A[train_index], A[test_index]
            nm_mse, nm_r2 = train_and_evaluate_nelder_mead(X_train, X_test, y_train, y_test, A_train, A_test, initial_params)
            scores['nelder_mead']['mse'].append(nm_mse)
            scores['nelder_mead']['r2'].append(nm_r2)
        else:
            # Linear model
            X_train_scaled, X_test_scaled = scale_features(X_train, X_test)
            linear_mse, linear_r2, linear_model = train_and_evaluate_regression(LinearRegression(), X_train_scaled, X_test_scaled, y_train, y_test)
            scores['linear']['mse'].append(linear_mse)
            scores['linear']['r2'].append(linear_r2)
            trained_models['linear'] = linear_model

            # XGBoost model
            xgb_mse, xgb_r2, xgb_model = train_and_evaluate_regression(XGBRegressor(n_estimators=100, learning_rate=0.1, random_state=42), X_train_scaled, X_test_scaled, y_train, y_test)
            scores['xgboost']['mse'].append(xgb_mse)
            scores['xgboost']['r2'].append(xgb_r2)
            trained_models['xgboost'] = xgb_model

            # Random Forest model (add these lines)
            rf_mse, rf_r2, rf_model = train_and_evaluate_regression(RandomForestRegressor(n_estimators=100, random_state=42), X_train_scaled, X_test_scaled, y_train, y_test)
            scores['random_forest']['mse'].append(rf_mse)
            scores['random_forest']['r2'].append(rf_r2)
            trained_models['random_forest'] = rf_model

    return calculate_mean_std_scores(scores), trained_models

def calculate_mean_std_scores(scores: Dict[str, Dict[str, list]]) -> Dict[str, float]:
    """
    Calculate mean and standard deviation of scores.

    Args:
        scores (Dict[str, Dict[str, list]]): Dictionary of scores.

    Returns:
        Dict[str, float]: Dictionary of mean and standard deviation of scores.
    """
    mean_std_scores = {}
    for model, metrics in scores.items():
        mean_std_scores[f'{model}_mse_mean'] = np.mean(metrics['mse']) if metrics['mse'] else np.nan
        mean_std_scores[f'{model}_mse_std'] = np.std(metrics['mse']) if metrics['mse'] else np.nan
        mean_std_scores[f'{model}_r2_mean'] = np.mean(metrics['r2']) if metrics['r2'] else np.nan
        mean_std_scores[f'{model}_r2_std'] = np.std(metrics['r2']) if metrics['r2'] else np.nan
    return mean_std_scores

def calculate_permutation_importance(model, X: pd.DataFrame, y: np.ndarray):
    """
    Calculate and print permutation importance for a given model.
    
    Args:
        model: Trained model
        X (pd.DataFrame): Feature matrix
        y (np.ndarray): Target values
    """
    perm_importance = permutation_importance(model, X, y, n_repeats=10, random_state=42)
    feature_importance = pd.DataFrame({
        'feature': X.columns,
        'importance': perm_importance.importances_mean
    }).sort_values('importance', ascending=False)
    
    print("\nPermutation Importance:")
    print(feature_importance)



def plot_learning_curves(X: pd.DataFrame, y: np.ndarray, name, model):
    """
    Plot learning curves for a given model.
    
    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        model: The model to train.
        initial_params (np.ndarray): Initial parameters for Nelder-Mead.
        A (np.ndarray): Asymptote values.
    """
    train_sizes = np.linspace(0.1, 0.9, 10)
    train_errors = []
    test_errors = []
    X, y = shuffle(X, y, random_state=42)
    for train_size in train_sizes:
        n_train = int(train_size * len(X))
        X_train, X_test = X.iloc[:n_train], X.iloc[n_train:]
        y_train, y_test = y[:n_train], y[n_train:]
        # X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=train_size, random_state=42)

        X_train_scaled, X_test_scaled = scale_features(X_train, X_test)
        train_mse, _, _ = train_and_evaluate_regression(model, X_train_scaled, X_train_scaled, y_train, y_train)
        test_mse, _, _ = train_and_evaluate_regression(model, X_train_scaled, X_test_scaled, y_train, y_test)

        train_errors.append(train_mse)
        test_errors.append(test_mse)

    plt.figure()
    plt.plot(train_sizes, train_errors, label='Training Error')
    plt.plot(train_sizes, test_errors, label='Validation Error')
    plt.xlabel('Training Set Size')
    plt.ylabel('Mean Squared Error')
    plt.title(f'{name} Learning Curves')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    return plt



def main():
    pars = [
        "nearest_mature", "age", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
        "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
        "cwd"
    ]

    # X, y, initial_params, _ = load_and_preprocess_data("./3_python_scripts/non_aggregated_100k_all.csv", pars)
    X, y, initial_params, _ = load_and_preprocess_data("./3_python_scripts/non_aggregated_5000.csv", pars)

    # results, trained_models = cross_validate(X, y, initial_params, run_nelder_mead=False)
    
    # # Print results for other models
    # print("\nResults for other models:")
    # for model in ['linear', 'xgboost', 'random_forest']:
    #     print(f"\n{model.capitalize()} Model Results:")
    #     print(f"Cross-validation MSE: {results[f'{model}_mse_mean']:.4f} (+/- {results[f'{model}_mse_std']:.4f})")
    #     print(f"Cross-validation R^2: {results[f'{model}_r2_mean']:.4f} (+/- {results[f'{model}_r2_std']:.4f})")
    #     # Calculate and print permutation importance for each model (except Nelder-Mead)
    #     print(f"\nPermutation Importance for {model.capitalize()} Model:")
    #     print(calculate_permutation_importance(trained_models[model], X, y))
    
    plot_learning_curves(X, y, 'xgboost', XGBRegressor(n_estimators=100, learning_rate=0.1, random_state=42))
    plot_learning_curves(X, y, 'linear', LinearRegression())
    plot_learning_curves(X, y,'random forest', RandomForestRegressor(n_estimators=100, random_state=42))
    plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # nelder_mead_pars = [par for par in pars if par != "nearest_mature"]

    # X, y, initial_params, A = load_and_preprocess_data("./3_python_scripts/non_aggregated_100k_all.csv", nelder_mead_pars)

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

def test():


    pass

if __name__ == "__main__":
    main()