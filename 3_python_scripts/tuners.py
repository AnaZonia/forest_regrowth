# scripts/hyperparameter_tuning.py
import numpy as np
import ray
from ray import tune
from ray.tune.schedulers import ASHAScheduler
from skopt import gp_minimize
from skopt.space import Real
from skopt.utils import use_named_args
from scipy.optimize import dual_annealing
from deap import base, creator, tools, algorithms
import random
from itertools import product
from tqdm import tqdm

def optimize_with_grid_search(fun, params):
    # Define the grid search space
    param_grid = {
        'B0': np.linspace(100, 160, 3),
        'theta': np.linspace(0.1, 10, 2),
        'coeff_0': np.linspace(0, 10, 2),  # age
        **{f'coeff_{i}': np.linspace(-10, 10, 2) for i in range(1, len(params) - 1)},
        f'coeff_{len(params)-1}': np.linspace(0, 10, 2)  # cwd
    }

    # Generate all combinations of parameters
    param_combinations = list(product(*param_grid.values()))

    best_mse = float('inf')
    best_params = None
   
    pbar = tqdm(total=len(param_combinations), desc="Grid Search Progress")
    
    # Iterate through all parameter combinations
    for i, params in enumerate(param_combinations, 1):
        all_params = np.array(params)
        mse = fun(params=all_params)

        if mse < best_mse:
            best_mse = mse
            best_params = all_params
            print(f"New best MSE: {best_mse} with params: {best_params}")

            # Update progress bar
        pbar.update(1)

        # Print progress every 5% of total combinations
        if i % max(1, len(param_combinations) // 20) == 0:
            print(f"\nTested {i}/{len(param_combinations)} combinations ({i/len(param_combinations):.1%})")

    pbar.close()

    # Convert best_params to the same format as other optimizers
    best_config = {
        'B0': best_params[0],
        'theta': best_params[1],
        **{f'coeff_{i}': best_params[i] for i in range(len(params)-1)}
    }

    return best_config

def optimize_with_ray_tune(fun, params):
    ray.init()

    def train_model(config):
        all_params = np.array([
            config["B0"],
            config["theta"],
            *[config[f"coeff_{i}"] for i in range(len(params))]
        ])
        
        # Objective function with Ray Tune's parameters
        mse = fun(params=all_params)
        ray.train.report({"mse": mse})
    
    config = {
        "B0": tune.uniform(60, 140),
        "theta": tune.uniform(0.1, 10),
        "coeff_0": tune.uniform(0, 10),  # age
        **{f"coeff_{i}": tune.uniform(-10, 10) for i in range(1, len(params)-1)},
        f"coeff_{len(params)-1}": tune.uniform(0, 10)  # cwd
    }


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
        num_samples=10,  # number of different parameter combinations to try
        verbose=0
    )
    
    best_config = analysis.get_best_config(metric="mse", mode="min")
    ray.shutdown()
    return best_config


def optimize_with_skopt(fun, params):
    # Define the search space
    space = [
        Real(60, 140, name='B0'),
        Real(0.1, 10, name='theta'),
        Real(0, 10, name='coeff_0'),  # age
        *[Real(-10, 10, name=f'coeff_{i}') for i in range(1, len(params) - 1)],
        Real(0, 10, name=f'coeff_{len(params)-1}')  # cwd
    ]

    # Define the objective function for skopt
    @use_named_args(space)
    def objective(**config):
        all_params = np.array([
            config['B0'],
            config['theta'],
            *[config[f'coeff_{i}'] for i in range(len(params))]
        ])
        
        return fun(params=all_params)

    # Run the optimization
    result = gp_minimize(
        objective,
        space,
        n_calls=50,  # number of evaluations of the objective function
        n_random_starts=10,  # number of random initial points
        random_state=42
    )

    # Extract the best parameters
    best_config = {dim.name: value for dim, value in zip(space, result.x)}
    return best_config





def optimize_with_simulated_annealing(fun, params):
    def objective(x):
        return fun(params=x)

    bounds = [(60, 140), (0.1, 10)] + [(0, 10)] + [(-10, 10)] * (len(params) - 2) + [(0, 10)]
    
    result = dual_annealing(objective, bounds, maxiter=1000, initial_temp=5230, restart_temp_ratio=2e-5, visit=2.62, accept=-5.0, maxfun=1e7, seed=42)
    
    best_params = {
        'B0': result.x[0],
        'theta': result.x[1],
        **{f'coeff_{i}': result.x[i+2] for i in range(len(params))}
    }
    return best_params




def optimize_with_genetic_algorithm(fun, params):
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    toolbox.register("attr_float", random.uniform, -10, 10)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, n=len(params)+2)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    def evalOneMax(individual):
        return fun(params=np.array(individual)),

    toolbox.register("evaluate", evalOneMax)
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=3)

    population = toolbox.population(n=100)
    ngen = 50
    
    algorithms.eaSimple(population, toolbox, cxpb=0.5, mutpb=0.2, ngen=ngen, verbose=False)

    best_ind = tools.selBest(population, 1)[0]
    best_params = {
        'B0': best_ind[0],
        'theta': best_ind[1],
        **{f'coeff_{i}': best_ind[i+2] for i in range(len(params))}
    }
    return best_params
