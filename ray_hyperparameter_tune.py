space = {
    "B0": tune.uniform(0, 100),
    "A": tune.uniform(0, 500),
    "theta": tune.uniform(0, 1),
    "sd": tune.uniform(0, 1),
}

space.update({key: tune.uniform(0, 100) for key in pars})

# Initialize Ray
ray.init()


# Define the objective function for Ray Tune
def objective(config):
    return negative_log_likelihood(config, data)


# def run_model(path):
#     data = import_data(path)  # Replace "path_to_your_data.csv" with the actual path to your CSV file


def short_trial_dirname_creator(trial):
    # Example: Use a shorter, hashed version of the trial name
    trial_name = str(trial.trial_id)
    return os.path.join("ray_results", trial_name)


# Assuming your objective function and config are defined
analysis = tune.run(
    objective,
    config=space,
    num_samples=100,
    metric="loss",
    mode="min",
    trial_dirname_creator=short_trial_dirname_creator,
)

# Get the best parameters
best_trial = analysis.get_best_trial("loss", "min", "last")
best_params = best_trial.config

# print("Best parameters:", best_params)
# print("Best loss:", best_trial.last_result["loss"])
