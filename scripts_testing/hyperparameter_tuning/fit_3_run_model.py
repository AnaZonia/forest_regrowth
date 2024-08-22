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
from sklearn.model_selection import train_test_split# Split the data into training and testing sets


# def main():


# # Define climatic parameters that change yearly
# biomes = ["amaz", "atla", "both"]

# # Define land-use history intervals to import four dataframes
# intervals = ["5yr", "10yr", "15yr", "all"]

datafile = pd.read_csv("./data/15yr_amaz.csv")

for col in datafile.columns:
    print(col)
    
    
# Using Skicit-learn to split data into training and testing sets
train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.25, random_state = 42)


# if __name__ == "__main__":


#     main()
