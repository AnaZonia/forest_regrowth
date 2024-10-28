import itertools
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from scipy.interpolate import BSpline

# Growth curve function with safe exponentiation to avoid overflow
def growth_curve(X, A, age, B0, betas, theta):
    k = np.dot(X, betas) * age
    k = np.clip(k, 1e-10, 7)  # Keep k within a reasonable range
    result = B0 + (A - B0) * (1 - np.exp(-k)) ** theta
    return result

# Wrapper for curve_fit using global bounds and initial_params
def fit_curve_fit(X, A, age, y, initial_params, bounds):
    # Fit the curve
    params, _ = curve_fit(
        lambda X_flat, B0, *betas_theta: growth_curve(
            X_flat.reshape(-1, X.shape[1]), A, age, B0, betas_theta[:-1], betas_theta[-1]
        ),
        X.ravel(),
        y,
        p0=initial_params,
        bounds=bounds
    )
    return params

# Wrapper for minimize using Nelder-Mead with -inf penalty for bounds
def fit_minimize(X, A, age, y, initial_params, bounds, method):

    # Loss function with -inf return for out-of-bounds values
    def loss(params):
        B0, *betas, theta = params

        # Check if parameters are out of bounds
        if not (bounds[0][0] <= B0 <= bounds[1][0]):
            return float('inf')
        for i, beta in enumerate(betas):
            if not (bounds[0][i + 1] <= beta <= bounds[1][i + 1]):
                return float('inf')
        if not (bounds[0][-1] <= theta <= bounds[1][-1]):
            return float('inf')

        # Calculate the sum of squared errors if within bounds
        y_pred = growth_curve(X, A, age, B0, betas, theta)
        result = np.sum((y - y_pred) ** 2)
        return result

    # Call minimize using Nelder-Mead
    result = minimize(loss, initial_params, method=method, tol=1e-1)
    return result.x

# Cross-validation with option for optimizer choice
def cross_validate(df, predictors, optimizer="curve_fit", initial_params=None):
    X = df[predictors].values
    y = df['agbd'].values
    A = df['nearest_mature_biomass'].values
    age = df['age'].values

    if initial_params is None:
        initial_params = [100] + [0] * len(predictors) + [1.0]
    
    bounds = ([0] + [-5] * len(predictors) + [0], 
              [150] + [5] * len(predictors) + [10])

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    r_squared_values = []

    for train_idx, test_idx in kf.split(X):
        scaler = MinMaxScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_test = scaler.transform(X[test_idx])
        y_train, y_test = y[train_idx], y[test_idx]
        A_train, A_test = A[train_idx], A[test_idx]
        age_train, age_test = age[train_idx], age[test_idx]

        # Choose the optimizer
        if optimizer == "curve_fit":
            params = fit_curve_fit(X_train, A_train, age_train, y_train, initial_params, bounds)
        else:
            params = fit_minimize(X_train, A_train, age_train, y_train, initial_params, bounds, optimizer)

        # Predict on the test set
        B0, *betas, theta = params
        y_pred = growth_curve(X_test, A_test, age_test, B0, betas, theta)
        print(np.dot(X_test, betas) * age_test)

        # Calculate R-squared
        ss_res = np.sum((y_test - y_pred) ** 2)
        ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        r_squared_values.append(r_squared)

    return np.mean(r_squared_values)


def run_model_iteration(data, param, best, base_row, all_pars_iter, categorical):
    if param in categorical:
        inipar = {**best['par'], **{key: all_pars_iter[key] for key in all_pars_iter.keys() if param in key}}
    else:
        inipar = {**best['par'], param: all_pars_iter[param]}
    
    model = minimize(lambda params: np.sum((data - params) ** 2), list(inipar.values()), method='Nelder-Mead')

    iter_row = base_row.copy()
    iter_row.update(dict(zip(inipar.keys(), model.x)))
    iter_row['likelihood'] = model.fun
    
    return iter_row

def iterative_model_selection(data, data_pars_iter, all_pars_iter, basic_pars_iter, categorical):
    remaining = list(range(len(data_pars_iter)))
    taken = len(remaining) + 1  # Ensure out-of-range index in the first iteration

    best = {'AIC': float('inf'), 'par': {key: all_pars_iter[key] for key in basic_pars_iter}}
    val = len(basic_pars_iter)

    base_row = {key: np.nan for key in all_pars_iter.keys()}
    base_row['likelihood'] = 0

    should_continue = True
    ideal_par_combination = []

    for i in range(len(data_pars_iter)):
        if not should_continue:
            break
        
        # Parallel processing for optimizing remaining parameters
        optim_remaining_pars = Parallel(n_jobs=-1)(delayed(run_model_iteration)(
            data, data_pars_iter[j], best, base_row, all_pars_iter, categorical) for j in remaining if j not in taken)

        # Create a DataFrame for results
        iter_df = pd.DataFrame(optim_remaining_pars)
        best_model_idx = iter_df['likelihood'].idxmin()
        best_model_AIC = 2 * iter_df['likelihood'].iloc[best_model_idx] + 2 * (i + val + 1)

        print(f"Iteration: {i + 1}, Parameters included: {i + 1}")

        # Update best model if AIC improves
        if best['AIC'] == float('inf') or best_model_AIC < best['AIC']:
            best['AIC'] = best_model_AIC
            best['par'] = iter_df.loc[best_model_idx, all_pars_iter.keys()].dropna().to_dict()
            taken = [idx for idx, name in enumerate(data_pars_iter) if any(name in k for k in best['par'].keys())]
        else:
            print("No improvement. Exiting loop.")
            should_continue = False

    ideal_par_combination.append(best['par'])
    return ideal_par_combination


# Cross-validation with option for optimizer choice, applied to each cluster
def cross_validate_cluster(df, predictors, optimizer="curve_fit", initial_params=None):
    # Apply PCA for dimensionality reduction and clustering to assign clusters
    X = df[predictors].values
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=15)  # Adjust based on explained variance
    X_pca = pca.fit_transform(X_scaled)
    kmeans = KMeans(n_clusters=15, init='k-means++', n_init=20, random_state=42)
    df['cluster'] = kmeans.fit_predict(X_pca)
    
    # Initialize a dictionary to store R-squared values per cluster
    cluster_r_squared = {}
    
    # Loop through each unique cluster and perform cross-validation
    for cluster_id in df['cluster'].unique():
        print(cluster_id)
        # Filter data for the current cluster
        cluster_data = df[df['cluster'] == cluster_id]
        
        # Define the predictors and response for the cluster
        X_cluster = cluster_data[predictors].values
        y_cluster = cluster_data['agbd'].values
        A_cluster = cluster_data['nearest_mature_biomass'].values
        age_cluster = cluster_data['age'].values

        # Set initial parameters and bounds if not provided
        if initial_params is None:
            initial_params = [100] + [0] * len(predictors) + [1.0]

        bounds = ([0] + [-5] * len(predictors) + [0], 
                  [150] + [5] * len(predictors) + [10])

        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        r_squared_values = []

        # Cross-validation loop within the cluster
        for train_idx, test_idx in kf.split(X_cluster):
            X_train = scaler.fit_transform(X_cluster[train_idx])
            X_test = scaler.transform(X_cluster[test_idx])
            y_train, y_test = y_cluster[train_idx], y_cluster[test_idx]
            A_train, A_test = A_cluster[train_idx], A_cluster[test_idx]
            age_train, age_test = age_cluster[train_idx], age_cluster[test_idx]

            # Choose the optimizer
            if optimizer == "curve_fit":
                params = fit_curve_fit(X_train, A_train, age_train, y_train, initial_params, bounds)
            else:
                params = fit_minimize(X_train, A_train, age_train, y_train, initial_params, bounds, optimizer)

            # Predict on the test set
            B0, *betas, theta = params
            y_pred = growth_curve(X_test, A_test, age_test, B0, betas, theta)

            # Calculate R-squared
            ss_res = np.sum((y_test - y_pred) ** 2)
            ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            print(r_squared)
            r_squared_values.append(r_squared)

        # Store the mean R-squared for the current cluster
        cluster_r_squared[cluster_id] = np.mean(r_squared_values)
        print(f"Cluster {cluster_id} - Mean R-squared: {cluster_r_squared[cluster_id]:.4f}")

    return cluster_r_squared

# Load the data
df = pd.read_csv("./0_data/non_aggregated_all.csv")
df = pd.get_dummies(df, columns=['topography', 'ecoreg', 'indig', 'protec', 'last_LU'], drop_first=True)
# df = df.sample(10000, random_state=42)
# predictors = [col for col in df.columns if col not in ['age', 'agbd', 'nearest_mature_biomass']]
predictors = ["nitro", "num_fires_before_regrowth", "sur_cover"]

# mean_r_squared_curve_fit = cross_validate(df, predictors, optimizer="curve_fit")
# print(f"Mean cross-validated R-squared using curve_fit: {mean_r_squared_curve_fit}")
# mean_r_squared_nelder_mead = cross_validate(df, predictors, optimizer="Nelder-Mead")
# print(f"Mean cross-validated R-squared using Nelder-Mead: {mean_r_squared_nelder_mead}")
mean_r_squared_LBFGSB = cross_validate(df, predictors, optimizer="L-BFGS-B")
print(f"Mean cross-validated R-squared using L-BFGS-B: {mean_r_squared_LBFGSB}")


# best_params_lbfgsb, best_r_squared_lbfgsb = grid_search_initial_params(df, predictors, optimizer="L-BFGS-B")
# print(f"Best initial parameters and R-squared using L-BFGS-B: {best_params_lbfgsb}, R^2 = {best_r_squared_lbfgsb}")


# # Run cross-validation for each cluster
# cluster_r_squared = cross_validate_cluster(df, predictors, optimizer="Nelder-Mead")

# # Print R-squared values per cluster
# for cluster_id, r2 in cluster_r_squared.items():
#     print(f"Cluster {cluster_id} - Mean R-squared: {r2:.4f}")











# Grid search over different starting values for B0, k, and theta
def grid_search_initial_params(df, predictors, optimizer="curve_fit"):
    # Define ranges for initial parameter guesses
    B0_range = [50, 100, 150]
    beta_range = [-5, 0, 5]
    theta_range = [0.5, 1.0, 1.5]

    best_r_squared = -np.inf
    best_params = None

    # Loop over each combination of initial parameter guesses
    for B0, beta, theta in itertools.product(B0_range, beta_range, theta_range):
        n_predictors = df[predictors].values.shape[1]
        initial_params = [B0] + [beta] * n_predictors + [theta]

        # Run cross-validation with the current initial parameters
        mean_r_squared = cross_validate(df, predictors, optimizer, initial_params)
        
        # Check if this combination yields a better R-squared
        if mean_r_squared > best_r_squared:
            best_r_squared = mean_r_squared
            best_params = (B0, beta, theta)
            print(f"New best R-squared: {best_r_squared} with params: {best_params}")
        else:
            print("no improvement")

    return best_params, best_r_squared




# # Function to calculate AIC given true values, predictions, and the number of parameters
# def calculate_aic(y_true, y_pred, num_params):
#     """Calculate AIC given true values, predictions, and the number of parameters."""
#     residual_sum_of_squares = np.sum((y_true - y_pred) ** 2)
#     n = len(y_true)
#     aic = n * np.log(residual_sum_of_squares / n) + 2 * num_params
#     return aic


# # Stepwise beta selection without cross-validation, using the whole dataset
# def stepwise_beta_selection(df, predictors, initial_params):
#     global bounds  # Assuming bounds are defined globally

#     X = df[predictors].values
#     y = df['agbd'].values
#     A = df['nearest_mature_biomass'].values
#     age = df['age'].values

#     n_predictors = X.shape[1]
#     all_betas = list(range(n_predictors))  # Indices of all possible beta parameters
#     selected_betas = []  # Start with no beta parameters selected
#     best_aic = np.inf
#     best_params = initial_params[:1] + [0] * n_predictors + [initial_params[-1]]  # Initialize with B0 and theta

#     while True:
#         best_beta = None
#         improvement = False

#         # Test each unselected beta parameter individually
#         for beta_idx in all_betas:
#             if beta_idx in selected_betas:
#                 continue  # Skip already selected betas

#             # Temporarily include the current beta in `selected_betas`
#             temp_selected_betas = selected_betas + [beta_idx]

#             # Construct initial params, setting unselected betas to zero
#             temp_params = [initial_params[0]]  # Start with B0
#             temp_params += [initial_params[1] if i in temp_selected_betas else 0 for i in range(n_predictors)]
#             temp_params += [initial_params[-1]]  # End with theta
#             bounds = ([0] + [-5] * len(temp_selected_betas) + [0], 
#               [150] + [5] * len(temp_selected_betas) + [10])
            
#             # Fit the model using the whole dataset with the current set of betas
#             fitted_params = fit_minimize(X, A, age, y, temp_params, bounds, "Nelder-Mead")

#             # Predict using this configuration to calculate AIC
#             B0, *betas, theta = fitted_params
#             y_pred = growth_curve(X, A, age, B0, betas, theta)
#             aic = calculate_aic(y, y_pred, num_params=len(temp_selected_betas) + 2)

#             # Check if this configuration yields a better AIC
#             if aic < best_aic:
#                 best_aic = aic
#                 best_beta = beta_idx
#                 best_params = fitted_params
#                 improvement = True

#         # If we found a new best beta, add it to selected_betas
#         if improvement:
#             selected_betas.append(best_beta)
#             print(f"Selected beta {best_beta} with AIC: {best_aic}")
#         else:
#             break  # No improvement, exit the loop

#     # Return the best parameters and AIC
#     return best_params, best_aic


# # Example usage:
# initial_params = [100] + [0] * len(predictors) + [1.0]
# best_params, best_aic = stepwise_beta_selection(df, predictors, initial_params)
# print(f"Best parameters: {best_params}")
# print(f"Best AIC: {best_aic}")
