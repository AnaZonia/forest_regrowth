


# Define a function to unpack parameters from a list
def unpack_parameters(params):
    return {
        "B0": params[0],
        "A": params[1],
        "theta": params[2],
        "sd": params[3],
        "age": params[4]
    }

# write the growth curve with yearly climatic data and permanent non-climatic data
def growth_curve(pars, data):
    r"""This function defines the growth function and parameter dictionary"""

    return (
        pars["B0"]
        + pars["A"] * (1 - np.exp(-pars["age"] * data["age"])) ** pars["theta"]
    )

# Modify the neg_log_likelihood function to unpack parameters
def neg_log_likelihood(params, data):
    params = unpack_parameters(params)
    likelihood = norm.logpdf(
        data["agbd"] - growth_curve(params, data), loc=0, scale=params["sd"]
    )
    return -np.sum(likelihood)

def fit_scipy(data):
    # Initial guesses
    initial_guess = [40, 80, 1.5, 1, 1]  # B0, A, theta, sd, age
    # Define bounds
    bounds = [(0, 100), (50, 400), (0, None), (0, None), (0, None)]  # B0, A, theta, sd
    # Call minimize with bounds
    result = minimize(neg_log_likelihood, initial_guess, args=(data,), bounds=bounds)

    return result.x


def objective(config, data):
    params = unpack_parameters([config["B0"], config["A"], config["theta"], config["sd"], config["age"]])
    return neg_log_likelihood(params, data)


def fit_ray(data):
    # Define the search space
    space = {
        "B0": tune.uniform(40, 100),
        "A": tune.uniform(100, 200),
        "theta": tune.uniform(0, 1),
        "sd": tune.uniform(3, 8),
        "age": tune.uniform(0, 2)
    }

    # Initialize Ray
    ray.init()

    # Run the optimizer
    analysis = tune.run(
        objective,
        data = data,
        config = space,
        num_samples = 100,
        metric = "loss",
        mode = "min",
    )

    # Get the best parameters
    best_trial = analysis.get_best_trial("loss", "min", "last")
    best_params = best_trial.config

    return best_params