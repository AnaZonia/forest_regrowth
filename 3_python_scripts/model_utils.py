"""
model_utils.py

This module contains utility functions and classes for performing regression analysis and optimization.
It includes functions for plotting learning curves, performing cross-validation, and optimizing models
using Nelder-Mead and Stan MCMC methods.

Functions:
- plot_learning_curves: Plot learning curves for a given model pipeline.
- regression_cv: Perform cross-validation for regression models.
- nelder_mead_B0_theta: Objective function for Nelder-Mead optimization with B0 and theta parameters.
- nelder_mead_lag: Objective function for Nelder-Mead optimization with lag parameters.
- process_fold: Process a single fold for cross-validation using Nelder-Mead optimization.
- cross_validate_nelder_mead: Cross-validate a model using Nelder-Mead optimization.
- optimize_with_pystan_mcmc: Optimize model parameters using PyStan MCMC.

"""

import numpy as np
import matplotlib.pyplot as plt
from functools import partial
import pandas as pd
import multiprocessing
from numpy.lib.recfunctions import drop_fields

from scipy.stats import norm
from scipy.optimize import minimize
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, cross_validate, GridSearchCV, learning_curve
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance

import stan

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



def regression_cv(X, y, model, unseen_data, name, param_grid = None):
    """
    Perform cross-validation for regression models.

    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        model: The regression model to train.
        unseen_data (DataSet): Unseen data for evaluation.
        name (str): Name of the model for the plot title.
        param_grid (dict): Parameter grid for hyperparameter tuning.

    Returns:
        tuple: Mean R2, standard deviation of R2, unseen R2, permutation importance, and learning curve figure.
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
                            scoring = 'neg_mean_squared_error', refit = False)
        grid_search.fit(X, y)

        best_params = grid_search.best_params_
        model = model.set_params(**best_params)

    # After performing grid search, perform cross-validation
    cv_results = cross_validate(model, X, y, cv = kf, 
                        scoring = ['neg_mean_squared_error'],
                        return_estimator = True)
    mse_scores = -cv_results['test_neg_mean_squared_error']
    best_index = np.argmin(-cv_results['test_neg_mean_squared_error'])
    best_model = cv_results['estimator'][best_index]

    # Perform permutation importance
    _, test_indices = splits[best_index]
    X_test = X.iloc[test_indices]
    y_test = y[test_indices]
    perm_importance = permutation_importance(best_model, X_test, y_test, n_repeats = 10, random_state = 42)

    # Calculate R2 from MSE
    r2_scores = 1 - (mse_scores / np.var(y))
    mean_r2 = np.mean(r2_scores)
    std_r2 = np.std(r2_scores)

    # Calculate r2 of the best model in the unseen dataset
    y_pred = best_model.predict(unseen_data.X)
    unseen_r2 = r2_score(unseen_data.y, y_pred)

    fig = plot_learning_curves(X, y, best_model, name)

    return mean_r2, std_r2, unseen_r2, perm_importance, fig, best_model


def nelder_mead_lag(params, X, y, A, fold = 1, return_predictions = False):
    """
    Objective function for Nelder-Mead optimization with lag parameters.

    Args:
        params (list): List of parameters [m_base, sd_base, sd, theta, coeffs].
        X (np.ndarray): Feature matrix.
        y (np.ndarray): Target values.
        A (np.ndarray): Asymptote values.
        random_state (int): Random seed for reproducibility.
        return_predictions (bool): Flag to return predictions instead of negative log-likelihood.

    Returns:
        float or np.ndarray: Negative log-likelihood or predictions.
    """
    m_base, sd_base, sd, theta, k0 = params[:5]
    coeffs = params[5:]

    age = X.age
    X = drop_fields(X, 'age')

    np.random.seed(fold)
    re_base = np.random.normal(loc=0, scale=1, size=1000)
    log_normal_scaled_base = np.exp(re_base * sd_base + m_base)
    
    m_results = np.zeros((len(y), len(log_normal_scaled_base)))
    y_pred_all = np.zeros((len(y), len(log_normal_scaled_base)))

    X_array = np.array([X[col_name] for col_name in X.dtype.names])
    pars_coeffs = X_array * coeffs[:, np.newaxis]

    # Vectorized operations
    k = pars_coeffs.T[:, np.newaxis] * (age[:, np.newaxis] + log_normal_scaled_base) + k0
    
    y_pred_all = A[:, np.newaxis] * (1 - np.exp(-k))**theta
    m_results = norm.pdf(y_pred_all - y[:, np.newaxis], scale=sd)

    mean_m_results = np.mean(m_results, axis = 1)
    # Replace zero values with epsilon
    mean_m_results[mean_m_results < 1e-10] = 1e-10
    neg_log_likelihood = np.round(-np.sum(np.log(mean_m_results)), decimals = 4)
    print(neg_log_likelihood)

    if return_predictions:
        return np.mean(y_pred_all, axis = 1)
    else:
        return neg_log_likelihood
