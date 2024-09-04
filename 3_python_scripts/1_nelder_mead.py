# scripts/nelder_mead_fit.py
from scipy.optimize import minimize
import numpy as np
import pandas as pd


# write the growth curve with yearly climatic data and permanent non-climatic data
def growth_curve(params, data):

    return (
        pars["B0"]
        * (1 - np.exp(-pars["age"] * data["age"])) ** pars["theta"]
    )

# Modify the neg_log_likelihood function to unpack parameters
def objective_function(params, X, y):
    predictions = chapman_richards_model(params, X)
    return np.mean((y - predictions) ** 2)



def fit_nelder_mead(data):
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

