# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              
#                 Forest Regrowth Model Fitting and Comparison                 
#                                                                              
#                            Ana Avila - August 2024                           
#                                                                              
#     This script runs and compares the Chapman-Richards growth curve fit      
#     through Scipy (Nelder Mead) and Random Forest
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from ray import tune
# from ray.tune.suggest.bayesopt import BayesOptSearch


def main():
    region = "countrywide"
    
    # Define climatic parameters that change yearly
    climatic_pars = ["prec", "si"]

    # Define parameters that do not correspond to data, used for functional form
    non_data_pars = ["k0", "B0_exp", "B0", "theta"]
    
    # Define land-use history intervals to import four dataframes
    intervals = ["5yr", "10yr", "15yr", "all"]

    datafiles = [f"./data/{region}_{interval}.csv" for interval in intervals]
    dataframes = [import_climatic_data(file) for file in datafiles]


if __name__ == "__main__":


    main()
