# scripts/nelder_mead_fit.py
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
from scipy.optimize import minimize

def load_and_preprocess_data(filepath):
    
    df = pd.read_csv(filepath)
    # Assuming df is your DataFrame
    pars = [
        "age", "nearest_mature", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
        "lulc_sum_40", "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
        "cwd"
    ]

    # Select only the specified columns
    X = df[pars]    
    y = df['agbd']

    # Create the new DataFrame with one row
    # Fill the row with values
    initial_params = np.zeros(len(pars) + 2)
    initial_params[0] = y.mean()
    initial_params[1] = 1

    return X, y, initial_params

# Define the objective function (growth curve + MSE)
def objective_function(pars, X, y):

    B0, theta = pars[:2]
    # The remaining elements are the coefficients for calculating k
    coeffs = pars[2:]
    # Calculate k using vectorized operations
    # Exclude 'nearest_mature' from X for this calculation
    X_for_k = X.drop(columns=['nearest_mature'])
    k = np.dot(X_for_k.values, coeffs) + 1

    nearest_mature = X['nearest_mature'].values
    y_pred = B0 + (nearest_mature - B0) * (1 - np.exp(-k))**theta

    return mean_squared_error(y, y_pred)


# Optimize the parameters using Nelder-Mead
def optimize_parameters(X, y, initial_params):
    result = minimize(lambda params: objective_function(params, X, y), 
                      initial_params, 
                      method='Nelder-Mead')
    return result.x

# Perform 5-fold cross-validation
def cross_validate(optimized_params, X, y):
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    mse_scores = []
    
    for train_index, test_index in kf.split(X):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        
        # Optimize parameters on training data
        train_params = optimize_parameters(X_train, y_train, optimized_params)
        
        # Evaluate on test data
        mse = objective_function(train_params, X_test, y_test)
        mse_scores.append(mse)
    
    return np.mean(mse_scores), np.std(mse_scores)


def main():
    X, y, initial_params = load_and_preprocess_data("0_data/processed_atla_all.csv")
    
    # Perform cross-validation using the initial parameters
    mean_mse, std_mse = cross_validate(X, y, initial_params)

    print(f"Cross-validation MSE: {mean_mse:.4f} (+/- {std_mse:.4f})")

if __name__ == "__main__":
    main()