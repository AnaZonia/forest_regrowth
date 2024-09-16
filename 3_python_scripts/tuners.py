# scripts/hyperparameter_tuning.py
import numpy as np
import ray
from ray import tune
from ray.tune.schedulers import ASHAScheduler

from skopt import gp_minimize
from skopt.space import Real
from skopt.utils import use_named_args

import stan

def optimize_with_ray_tune(fun, pars):
    ray.init()

    def train_model(config):
        all_pars = np.array([
            config["B0"],
            config["theta"],
            *[config[f"coeff_{i}"] for i in range(len(pars))]
        ])
        
        # Objective function with Ray Tune's parameters
        mse = fun(pars=all_pars)
        
        ray.train.report({"mse": mse})
    
    config = {
        "B0": tune.uniform(60, 140),
        "theta": tune.uniform(0.1, 10)
    }

    # Add individual coefficient parameters
    config["coeff_0"] = tune.uniform(0, 10)# age
    for i in range(1,len(pars)-1):
        config[f"coeff_{i}"] = tune.uniform(-10, 10)
    config[f"coeff_{len(pars)-1}"] = tune.uniform(0, 10) # cwd


    scheduler = ASHAScheduler(
        metric="mse",
        mode="min",
        max_t=100,  # maximum number of training iterations
        grace_period=1,
        reduction_factor=2
    )
    
    analysis = tune.run(
        train_model,
        config=config,
        scheduler=scheduler,
        num_samples=10  # number of different parameter combinations to try
    )
    
    best_config = analysis.get_best_config(metric="mse", mode="min")
    ray.shutdown()
    return best_config

def optimize_with_skopt(fun, pars):
    # Define the search space
    space = [
        Real(60, 140, name='B0'),
        Real(0.1, 10, name='theta'),
        Real(0, 10, name='coeff_0'),  # age
    ]
    
    # Add individual coefficient parameters
    for i in range(1, len(pars) - 1):
        space.append(Real(-10, 10, name=f'coeff_{i}'))
    
    space.append(Real(0, 10, name=f'coeff_{len(pars)-1}'))  # cwd

    # Define the objective function for skopt
    @use_named_args(space)
    def objective(**params):
        all_pars = np.array([
            params['B0'],
            params['theta'],
            *[params[f'coeff_{i}'] for i in range(len(pars))]
        ])
        
        return fun(pars=all_pars)

    # Run the optimization
    result = gp_minimize(
        objective,
        space,
        n_calls=50,  # number of evaluations of the objective function
        n_random_starts=10,  # number of random initial points
        random_state=42
    )

    # Extract the best parameters
    best_params = {
        'B0': result.x[0],
        'theta': result.x[1],
    }
    for i in range(len(pars)):
        best_params[f'coeff_{i}'] = result.x[i+2]

    return best_params


def optimize_with_pystan_mcmc(fun, pars, X, y, A):
    # Define the Stan model
    stan_model = """
    data {
        int<lower=0> N;
        int<lower=0> P;
        matrix[N, P] X;
        vector[N] y;
        vector[N] A;
    }
    parameters {
        real<lower=0, upper=200> B0;
        real<lower=0, upper=10> theta;
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
    sm = pystan.StanModel(model_code=stan_model)

    # Fit the model
    fit = sm.sampling(data=stan_data, iter=2000, chains=4, warmup=1000, n_jobs=-1)

    # Extract the posterior means
    B0 = fit['B0'].mean()
    theta = fit['theta'].mean()
    coeffs = fit['coeffs'].mean(axis=0)

    # Combine parameters
    optimized_params = np.concatenate(([B0, theta], coeffs))

    return optimized_params

# Example usage:
# best_params = optimize_with_pystan_mcmc(objective_function, pars, X, y, A)
