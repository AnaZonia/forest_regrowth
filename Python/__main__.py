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
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold

from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

df = pd.read_csv("data/processed_amaz_10yr.csv")
df = df.drop(columns=[col for col in df.columns if col.startswith('prec_') or col.startswith('si_')])
# for col in df.columns:
#     print(col)
print(df.shape)
# Separate the response variable (agbd) and predictors

df['ecoreg'] = df['ecoreg'].astype('category')
df['biome'] = df['biome'].astype('category')
df['indig'] = df['indig'].astype('category')
df['protec'] = df['protec'].astype('category')
df['soil'] = df['soil'].astype('category')

X = df.drop(columns=['agbd', 'Unnamed: 0', 'biome', 'distance'])
# X = df[['cwd', 'mean_prec', 'mean_si', 'ecoreg', 'nearest_mature']]
y = df['agbd']                 # Response variable
print(df.dtypes)
# Fit the Random Forest Regressor
# Create the Random Forest Regressor
model = RandomForestRegressor(n_estimators=200, random_state=42)

# Create the KFold object
kf = KFold(n_splits=5, shuffle=True, random_state=42)

# Perform 5-fold cross-validation
cv_scores = cross_val_score(model, X, y, cv=kf, scoring='r2')

# Print the results
print("Cross-validation scores (R-squared):")
for i, score in enumerate(cv_scores):
    print(f"Fold {i+1}: {score:.4f}")

print(f"\nMean R-squared: {cv_scores.mean():.4f}")
print(f"Standard deviation of R-squared: {cv_scores.std():.4f}")

# Fit the model on the entire dataset to access feature importances
model.fit(X, y)

# Print feature importance scores
feature_importances = model.feature_importances_
print("\nFeature importance scores:")
for i, importance in enumerate(feature_importances):
    print(f"Feature {i+1} ({X.columns[i]}): {importance:.4f}")

# print('atlantic')
# # Calculate permutation importance
# perm_importance = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=42)

# # Sort the importance values
# sorted_importance_idx = perm_importance.importances_mean.argsort()
# for idx in sorted_importance_idx:
#     print(f"Feature: {X.columns[idx]}, Importance: {perm_importance.importances_mean[idx]}")


# # # Predict on the test set
# y_pred = model.predict(X_test)

# # # Evaluate the model
# # mse = mean_squared_error(y_test, y_pred)
# r2 = r2_score(y_test, y_pred)

# # print(f"Mean Squared Error: {mse}")
# print(f"R-squared: {r2}")

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
