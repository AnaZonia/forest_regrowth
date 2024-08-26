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
from sklearn.inspection import permutation_importance

df = pd.read_csv("data/15yr_amaz.csv")

df = df.drop(columns=[col for col in df.columns if col.startswith('prec_') or col.startswith('si_')])
# for col in df.columns:
#     print(col)
# print(df.shape)
X = df.drop(columns=[col for col in df.columns if 'ecoreg' in col])
# Separate the response variable (agbd) and predictors
X = X.drop(columns=['agbd', 'distance', 'last_LU_48'])  # Predictors (all columns except 'agbd')
y = df['agbd']                 # Response variable

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Fit the Random Forest Regressor
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)


# Calculate permutation importance
perm_importance = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=42)

# Sort the importance values
sorted_importance_idx = perm_importance.importances_mean.argsort()
for idx in sorted_importance_idx:
    print(f"Feature: {X.columns[idx]}, Importance: {perm_importance.importances_mean[idx]}")


# # Predict on the test set
y_pred = model.predict(X_test)

# # Evaluate the model
# mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

# print(f"Mean Squared Error: {mse}")
print(f"R-squared: {r2}")

# from xgboost import XGBRFRegressor
# from sklearn.model_selection import RepeatedKFold, cross_val_score
# from numpy import mean, std

# # define the model
# model = XGBRFRegressor(n_estimators=100, subsample=0.9, colsample_bynode=0.2)

# # define the model evaluation procedure
# cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)

# # evaluate the model and collect the scores
# # Change the scoring parameter to 'r2'
# r2_scores = cross_val_score(model, X, y, scoring='r2', cv=cv, n_jobs=-1)

# # report performance
# print('R-squared: %.3f (%.3f)' % (mean(r2_scores), std(r2_scores)))



# import pandas as pd
# from sklearn.model_selection import train_test_split
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.metrics import r2_score

# def iterative_feature_selection(X, y, random_state=42):
#     # Initialize variables
#     selected_predictors = []
#     best_r2 = -float('inf')
#     current_r2 = 0

#     # Create a list of all predictors
#     remaining_predictors = list(X.columns)

#     # Loop until R-squared decreases
#     while remaining_predictors:
#         r2_scores = []

#         # Try adding each remaining predictor one by one
#         for predictor in remaining_predictors:
#             predictors_to_test = selected_predictors + [predictor]
#             X_train, X_test, y_train, y_test = train_test_split(X[predictors_to_test], y, test_size=0.2, random_state=random_state)
#             model = RandomForestRegressor(n_estimators=100, random_state=random_state)
#             model.fit(X_train, y_train)
#             y_pred = model.predict(X_test)
#             r2 = r2_score(y_test, y_pred)
#             r2_scores.append((r2, predictor))

#         # Find the predictor with the highest R-squared
#         r2_scores.sort(reverse=True)
#         best_new_r2, best_predictor = r2_scores[0]

#         # Stop if R-squared decreases
#         if best_new_r2 > best_r2:
#             best_r2 = best_new_r2
#             selected_predictors.append(best_predictor)
#             remaining_predictors.remove(best_predictor)
#             print(f"Selected Predictors: {selected_predictors}, R-squared: {best_r2}")
#         else:
#             print(f"R-squared decreased. Stopping with predictors: {selected_predictors}")
#             break

#     return selected_predictors, best_r2

# # Usage example
# X = df.drop(columns=['agbd', 'distance', 'last_LU_48'])  # Predictors (all columns except 'agbd')
# y = df['agbd']                 # Response variable

# selected_predictors, best_r2 = iterative_feature_selection(X, y)
# print(f"Final selected predictors: {selected_predictors}")
# print(f"Best R-squared value: {best_r2}")
