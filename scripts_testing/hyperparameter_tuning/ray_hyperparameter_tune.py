import ray
from ray import tune

# Assuming the other functions and imports are defined as before


def objective(config):
    params = unpack_parameters(config)
    likelihood = norm.logpdf(
        data["agbd"] - growth_curve(params, data), loc=0, scale=params["sd"]
    )
    return -np.sum(likelihood)


def main():
    path = "./data/unified_data_15_years.csv"
    data = import_data(path)

    # Define the search space
    space = {
        "B0": tune.uniform(0, 100),
        "A": tune.uniform(0, 500),
        "theta": tune.uniform(0, 1),
        "sd": tune.uniform(0, 1),
    }

    # Initialize Ray
    ray.init()

    # Run the optimizer
    analysis = tune.run(
        objective,
        config=space,
        num_samples=100,
        metric="loss",
        mode="min",
    )

    # Get the best parameters
    best_trial = analysis.get_best_trial("loss", "min", "last")
    best_params = best_trial.config

    print(best_params)


if __name__ == "__main__":
    main()