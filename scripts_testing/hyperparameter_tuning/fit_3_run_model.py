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
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

# def main():

# # Define climatic parameters that change yearly
# biomes = ["amaz", "atla", "both"]

# # Define land-use history intervals to import four dataframes
# intervals = ["5yr", "10yr", "15yr", "all"]

df = pd.read_csv("data/all_amaz.csv")

df = df.drop(columns=[col for col in df.columns if col.startswith('prec_') or col.startswith('si_')])
# for col in df.columns:
#     print(col)
print(df.shape)

# Separate the response variable (agbd) and predictors
X = df.drop(columns=['agbd', 'distance', 'last_LU_48'])  # Predictors (all columns except 'agbd')
y = df['agbd']                 # Response variable

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Fit the Random Forest Regressor
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Predict on the test set
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error: {mse}")
print(f"R-squared: {r2}")