import numpy as np
import pandas as pd

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




