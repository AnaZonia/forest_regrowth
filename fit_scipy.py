r"""
Date 
    05/05/2024

Purpose
    This script runs. Hopefully.
    
Author
    Ana Catarina Avila
    McGill University
"""

import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("./forest_regrowth/")


def import_data(path: str):
    r"""This function:
    - Imports the dataframe
    - Removes unnecessary columns that will not be used in analysis
    - Converts categorical data to dummy variables"""

    # Read the CSV file
    data = pd.read_csv(path)

    # Drop unnecessary columns
    data = data.drop(columns=["system:index", ".geo", "biome"])

    # Convert 'soil' and 'ecoreg' to categorical data
    categorical = ["soil", "ecoreg"]
    data[categorical] = data[categorical].astype("category")

    # Create dummy variables for 'ecoreg' and 'soil'
    data = pd.get_dummies(data, columns=["ecoreg", "soil"])

    return data


data = import_data("./data/final_w_biomass.csv")

old_data = pd.read_csv("./data/amazon_df_sample_10million.csv")

# Define the years and climatic variables
years = np.arange(1985, 2020)

# define the climatic parameters - the ones that change yearly
climatic = ["prec", "si"]
# define the non-climatic parameters - the ones that are fixed throughout regrowth and
# that are used for fitting the model (excludes age and agbd)
non_climatic = [
    col for col in data.columns if not any(var in col for var in climatic + ["agbd"])
]

# define the dictionary that will hold the parameters and their values
pars = {
    "B0": 40,
    "A": 80,
    "theta": 1.5,
    "sd": 1,
}  # intercept, asymptote, shape term, standard deviation
pars_init = climatic + non_climatic
pars_init_values = [0.1] * len(pars_init)
new_elements = dict(zip(pars_init, pars_init_values))
pars.update(new_elements)


# write the growth curve with yearly climatic data and permanent non-climatic data
def growth_curve(pars, data, years):
    r"""This function defines the growth function and parameter dictionary"""

    # Calculate k
    k = np.zeros(len(data))
    for year in years:
        for clim_var in climatic:
            k += pars[clim_var] * data[f"{clim_var}_{year}"]
        for unique_var in non_climatic:
            k += pars[unique_var] * data[unique_var]

    return (
        pars["B0"]
        + pars["A"]
        * (1 - np.exp(-pars["last_fire"] * data["last_fire"])) ** pars["theta"]
    )


def negative_log_likelihood(pars, data):
    # Assuming a normal distribution for the data
    likelihood = norm.pdf(
        data["agbd"] - growth_curve(pars, data, years), loc=0, scale=pars["sd"]
    )
    return -np.sum(np.log(likelihood))


x0 = list(pars.values())


# Now you can call minimize with the modified growth_curve function
result = minimize(negative_log_likelihood, x0, args=(data,), method="Nelder-Mead")


# if __name__=="__main__":
#     path = './data/final_w_biomass.csv'
#     years = np.arange(1985, 2020)
#     # Initialize parameters
#     pars = {'B0': 40, 'A': 80, 'age': 0.05, 'theta': 1.5}

#     # Import the data
#     run_model()
