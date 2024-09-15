# scripts/hyperparameter_tuning.py

import numpy as np

def optimize_with_ray_tune(fun):

    def train_model(config):
        all_pars = np.array([
            config["B0"],
            config["theta"],
            *[config[f"coeff_{i}"] for i in range(len(model_pars) - 2)]
        ])
        
        # Objective function with Ray Tune's parameters
        mse = fun(pars=all_pars)
        
        tune.report(mse=mse)
    
    config = {
        "B0": tune.uniform(0, 10),
        "theta": tune.uniform(0.1, 10)
    }

        # Add individual coefficient parameters
    for i in range(len(model_pars) - 2):
        config[f"coeff_{i}"] = tune.uniform(-1, 1)
    
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
    return best_config