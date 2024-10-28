import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
import ray
from ray import tune

# Growth curve function, k is multiplied by age
def growth_curve(X, A, age, B0, betas, theta):
    k = np.dot(X, betas) * age
    k = np.clip(k, 1e-10, 7)  # Constrain k to [1e-10, 7]
    return B0 + (A - B0) * (1 - np.exp(-k)) ** theta

# Loss function
def loss(params, X, y, A, age):
    B0, *betas, theta = params
    y_pred = growth_curve(X, A, age, B0, np.array(betas), theta)
    return np.sum((y - y_pred) ** 2)

# Gradient approximation for SGD using finite differences
def compute_gradient(loss_fn, params, X, y, A, age, epsilon=1e-5):
    gradients = np.zeros_like(params)
    for i in range(len(params)):
        params_eps = np.array(params, dtype=float)
        params_eps[i] += epsilon
        loss1 = loss_fn(params_eps, X, y, A, age)
        
        params_eps[i] -= 2 * epsilon
        loss2 = loss_fn(params_eps, X, y, A, age)
        
        gradients[i] = (loss1 - loss2) / (2 * epsilon)
    return gradients

# SGD optimizer
def stochastic_gradient_descent(loss_fn, initial_params, X, y, A, age, lr=0.001, n_iter=1000):
    params = np.array(initial_params, dtype=float)
    for _ in range(n_iter):
        gradients = compute_gradient(loss_fn, params, X, y, A, age)
        params -= lr * gradients
    return params

# Cross-validation function with hyperparameter tuning using Ray Tune
def cross_validate_with_tuning(df, predictors):
    def train_model(config):
        X = df[predictors].values
        y = df['agbd'].values
        A = df['nearest_mature_biomass'].values
        age = df['age'].values

        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        r_squared_values = []

        for train_idx, test_idx in kf.split(X):
            scaler = MinMaxScaler()
            X_train = scaler.fit_transform(X[train_idx])
            X_test = scaler.transform(X[test_idx])
            y_train, y_test = y[train_idx], y[test_idx]
            A_train, A_test = A[train_idx], A[test_idx]
            age_train, age_test = age[train_idx], age[test_idx]

            n_predictors = X_train.shape[1]
            initial_params = [config["B0"]] + [config["betas"]] * n_predictors + [config["theta"]]

            if config["optimizer"] == "nelder-mead":
                result = minimize(loss, initial_params, args=(X_train, y_train, A_train, age_train), method='Nelder-Mead')
                params = result.x
            elif config["optimizer"] == "sgd":
                params = stochastic_gradient_descent(loss, initial_params, X_train, y_train, A_train, age_train,
                                                     lr=config["lr"], n_iter=config["n_iter"])

            B0, *betas, theta = params
            y_pred = growth_curve(X_test, A_test, age_test, B0, np.array(betas), theta)

            ss_res = np.sum((y_test - y_pred) ** 2)
            ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            r_squared_values.append(r_squared)

        tune.report(mean_r_squared=np.mean(r_squared_values))

    search_space = {
        "B0": tune.grid_search([50, 100, 150]),
        "betas": tune.grid_search([0.1, 0.5, 1.0]),
        "theta": tune.grid_search([0.5, 1.0, 2.0]),
        "optimizer": tune.grid_search(["nelder-mead", "sgd"]),
        "lr": tune.choice([0.0001, 0.001, 0.01]),  # Learning rate for SGD
        "n_iter": tune.choice([500, 1000, 1500])  # Iterations for SGD
    }

    ray.init(ignore_reinit_error=True)
    analysis = tune.run(train_model, config=search_space, metric="mean_r_squared", mode="max", num_samples=1)
    ray.shutdown()

    best_config = analysis.get_best_config(metric="mean_r_squared", mode="max")
    best_score = analysis.best_result["mean_r_squared"]
    
    return best_config, best_score

# Load the data
df = pd.read_csv("./0_data/non_aggregated_all.csv")        
df = pd.get_dummies(df, columns=['topography', 'ecoreg', 'indig', 'protec', 'last_LU'], drop_first=True)
df = df.sample(10000, random_state=42)

# Example usage
predictors = [col for col in df.columns if col not in ['age', 'agbd', 'nearest_mature_biomass']]
best_config, best_score = cross_validate_with_tuning(df, predictors)
print(f"Best configuration: {best_config}")
print(f"Best mean R-squared: {best_score}")
