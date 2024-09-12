# main.py
from data_processing import load_and_preprocess_data
from nelder_mead import objective_function, cross_validate
from scipy.optimize import minimize

def main():
    X, y, initial_params = load_and_preprocess_data("0_data/processed_atla_all.csv")
    
    # Optimize using Nelder-Mead
    result = minimize(cross_validate, initial_params, method='Nelder-Mead')

    # Fit Nelder-Mead
    # nelder_mead_params = fit_nelder_mead(X, y)
    # print(f"Nelder-Mead fitted parameters: {nelder_mead_params}")

    # # Fit Random Forest
    # rf_model, cv_scores, feature_importances = fit_random_forest(X, y)
    # print(f"Random Forest CV scores: {cv_scores}")
    # print(f"Random Forest feature importances: {feature_importances}")

    # # Hyperparameter tuning
    # best_model, best_params, best_score = bayesian_optimization(X, y)
    # print(f"Best model parameters: {best_params}, Best R2 score: {best_score}")

    # # Evaluate the best model
    # mse, r2 = evaluate_model(best_model, X, y)
    # print(f"Evaluation - MSE: {mse}, R2: {r2}")

if __name__ == "__main__":
    main()





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
