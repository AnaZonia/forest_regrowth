import numpy as np
import matplotlib.pyplot as plt
from functools import partial
import pandas as pd

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

    return mean_r2, std_r2, unseen_r2, perm_importance, fig


def nelder_mead_B0_theta(params, X, y, A, return_predictions = False):
    B0, theta = params[:2]
    coeffs = params[2:]
    
    adjustment_value = -np.log(1 - (y.mean() / A.mean()))
    k = np.dot(X, coeffs) + adjustment_value
    k = np.where(k < 0, adjustment_value, k)

    y_pred = B0 + (A - B0) * (1 - np.exp(-k))**theta
    mse = np.mean((y - y_pred)**2)

    if return_predictions:
        return y_pred
    else:
        return mse

def nelder_mead_lag(params, X, y, A, random_state = 1, return_predictions=False):

    m_base, sd_base, sd, theta = params[:4]
    coeffs = params[4:]
    
    age = X['age']
    X = X.drop('age', axis = 1)

    np.random.seed(random_state)
    re_base = np.random.normal(0, 1, 1000)
    log_normal_scaled_base = np.exp((re_base + m_base) * sd_base)
    
    m_results = np.zeros((len(y), len(log_normal_scaled_base)))
    y_pred_all = np.zeros((len(y), len(log_normal_scaled_base)))

    for i, lag in enumerate(log_normal_scaled_base):
        k = np.dot(X, coeffs) * age + lag
        y_pred = A * (1 - np.exp(-k))**theta
        y_pred_all[:, i] = y_pred
        m_results[:, i] = norm.pdf(y_pred - y, scale = sd)

    mean_m_results = np.mean(m_results, axis = 1)
    # Replace zero values with epsilon
    mean_m_results[mean_m_results < 1e-10] = 1e-10
    neg_log_likelihood = -np.sum(np.log(mean_m_results))
    print(neg_log_likelihood)

    if return_predictions:
        return np.mean(y_pred_all, axis = 1)
    else:
        return neg_log_likelihood


def cross_validate_nelder_mead(X, y, A, params, unseen_data, name):
    fold_scores = []
    fit_params = []
    scaler = MinMaxScaler()
    kf = KFold(n_splits = 5, shuffle = True, random_state = 42)

    for fold, (train_index, test_index) in enumerate(kf.split(X)):
        print(fold)
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y[train_index], y[test_index]
        A_train, A_test = A[train_index], A[test_index]

        # Fit and transform X_train, keeping it as a DataFrame
        X_train = pd.DataFrame(
            scaler.fit_transform(X_train),
            columns=X_train.columns,
            index=X_train.index
        )

        # Transform X_test, keeping it as a DataFrame
        X_test = pd.DataFrame(
            scaler.transform(X_test),
            columns=X_test.columns,
            index=X_test.index
        )

        # Using fold as a random state for reproducibility
        result = minimize(nelder_mead_lag, params, args = (X_train, y_train, A_train, fold), \
                          method = 'Nelder-Mead')

        val_score = nelder_mead_lag(result.x, X_test, y_test, A_test, fold)
        fit_params.append(result.x)
        fold_scores.append(val_score)

    mean_score = np.mean(fold_scores)
    std_score = np.std(fold_scores)
    # Find the best model
    best_index = np.argmin(fold_scores)
    best_params = fit_params[best_index]

    # Calculate r2 of the best model in the unseen dataset
    unseen_X_scaled = scaler.fit_transform(unseen_data.X)
    y_pred = nelder_mead_lag(best_params, unseen_X_scaled, unseen_data.y, unseen_data.A, \
                         return_predictions = True)
    unseen_r2 = r2_score(unseen_data.y, y_pred)

    # model_pipeline = Pipeline([
    #     ('scaler', MinMaxScaler()),
    #     ('regressor', partial(nelder_mead_lag, A = A))
    # ])

    # fig = plot_learning_curves(X, y, model_pipeline, name)

    return mean_score, std_score, unseen_r2



def optimize_with_pystan_mcmc(fun, pars, X, y, A):
    # Define the Stan model
    stan_model = """
    data {
        int<lower = 0> N;
        int<lower = 0> P;
        matrix[N, P] X;
        vector[N] y;
        vector[N] A;
    }
    parameters {
        real<lower = 0, upper = 200> B0;
        real<lower = 0, upper = 10> theta;
        vector[P] coeffs;
    }
    model {
        vector[N] mu;
        vector[N] k;
        
        // Priors
        B0 ~ uniform(60, 140);
        theta ~ uniform(0.1, 10);
        coeffs ~ normal(0, 10);
        
        // Model
        k = X * coeffs + (-log1m(mean(y) / mean(A)));
        for (n in 1:N) {
            if (k[n] < 0) {
                k[n] = -log1m(mean(y) / mean(A));
            }
        }
        mu = B0 + (A - B0) .* (1 - exp(-k)) .^ theta;
        
        // Likelihood
        y ~ normal(mu, 1);  // Assuming unit variance for simplicity
    }
    """

    # Prepare data for Stan
    stan_data = {
        'N': X.shape[0],
        'P': X.shape[1],
        'X': X,
        'y': y,
        'A': A
    }

    # Compile the model
    sm = stan.StanModel(model_code = stan_model)

    # Fit the model
    fit = sm.sampling(data = stan_data, iter = 2000, chains = 4, warmup = 1000, n_jobs = -1)

    # Extract the posterior means
    B0 = fit['B0'].mean()
    theta = fit['theta'].mean()
    coeffs = fit['coeffs'].mean(axis = 0)

    # Combine parameters
    optimized_params = np.concatenate(([B0, theta], coeffs))

    return optimized_params

# Example usage:
# best_params = optimize_with_pystan_mcmc(objective_function, pars, X, y, A)
