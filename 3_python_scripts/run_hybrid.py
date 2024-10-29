import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize
from sklearn.cluster import KMeans
from sklearn.model_selection import KFold, train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.metrics import r2_score
# Global initial guess and bounds for parameters [B0, betas..., theta]
n_predictors = None  # Define this later based on predictors
initial_params = None
bounds = None

# Load the data
df = pd.read_csv("./0_data/non_aggregated_all.csv")
df = pd.get_dummies(df, columns=['topography', 'ecoreg', 'indig', 'protec', 'last_LU'], drop_first=True)
# df = df.sample(10000, random_state=42)

# Example usage for comparison
predictors = [col for col in df.columns if col not in ['age', 'agbd', 'nearest_mature_biomass']]

# Normalize predictor values
scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(df[predictors])

# Step 2: Apply PCA for dimensionality reduction
pca = PCA(n_components=15)  # Adjust based on explained variance
X_pca = pca.fit_transform(X_scaled)
kmeans = KMeans(n_clusters=15, init='k-means++', n_init=20, random_state=42)
kmeans.fit(X_pca)
df['cluster'] = kmeans.fit_predict(X_pca)

# Growth curve function with safe exponentiation to avoid overflow
def growth_curve_k(A, age, k, B0):
    return B0 + (A - B0) * (1 - np.exp(-k * age))

# Wrapper for minimize using Nelder-Mead with -inf penalty for bounds
def fit_minimize(A, age, y, initial_params, bounds, method):

    # Loss function with -inf return for out-of-bounds values
    def loss(params):
        B0, k = params
        # Check if parameters are out of bounds
        if not (bounds[0][0] <= B0 <= bounds[1][0]):
            return float('inf')
        # Check if parameters are out of bounds
        if not (bounds[0][1] <= k <= bounds[1][1]):
            return float('inf')

        # Calculate the sum of squared errors if within bounds
        y_pred = growth_curve_k(A, age, k, B0)
        return np.sum((y - y_pred) ** 2)

    # Call minimize using Nelder-Mead
    result = minimize(loss, initial_params, method=method)
    return result.x

# Initialize a list to store the 'k' values for each row
k_values = []

# Step 4: Fit the growth curve for each cluster and extract 'k' parameter
for cluster_id in df['cluster'].unique():
    # Filter data for the current cluster
    cluster_data = df[df['cluster'] == cluster_id]
    
    # Extract age and agbd columns for curve fitting
    age = cluster_data['age'].values
    A = cluster_data['nearest_mature_biomass'].values
    agbd = cluster_data['agbd'].values  # The target values for fitting

    # Set initial guesses and bounds for B0, k, and theta
    initial_guess = [100, 1]
    bounds = ([0, 0], [150, 7])

    # Fit the growth curve for the current cluster
    params = fit_minimize(A, age, agbd, initial_guess, bounds, method='Nelder-Mead')

    B0_value, k_value = params

    # Calculate predictions and R-squared for the cluster
    y_pred = growth_curve_k(A, age, k_value, B0_value)
    ss_res = np.sum((agbd - y_pred) ** 2)
    ss_tot = np.sum((agbd - np.mean(agbd)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    # Print R-squared for the cluster
    print(f"Cluster {cluster_id} - B0: {B0_value:.2f}, k: {k_value:.4f}, R-squared: {r_squared:.4f}")

    # Assign 'k' value to all rows in the cluster
    k_values.extend([k_value] * len(cluster_data))

# Add the 'k' values as a new column to the DataFrame
df['k'] = k_values
# Separate predictors and target
X_scaled = df[predictors].values  # Assuming predictors contains the columns for X
y = df['k'].values  # The target variable is now 'k'

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# Initialize and train the XGBoost regressor
xgb_regressor = xgb.XGBRegressor(objective='reg:squarederror', n_estimators=100, max_depth=5, learning_rate=0.1, random_state=42)
xgb_regressor.fit(X_train, y_train)

# Predict on the test set
y_pred = xgb_regressor.predict(X_test)

r2 = r2_score(y_test, y_pred)

print(f"R-squared (R2): {r2}")

df['k_predicted'] = xgb_regressor.predict(X_scaled)

def growth_curve(A, age, k, B0):
    result = B0 + (A - B0) * (1 - np.exp(-(k * age)))
    return result

# Wrapper for minimize using Nelder-Mead with -inf penalty for bounds
def fit_minimize(A, age, k, y, initial_params, bounds, method):

    # Loss function with -inf return for out-of-bounds values
    def loss(params):
        B0, k0 = params

        # Check if parameters are out of bounds
        if not (bounds[0][0] <= B0 <= bounds[1][0]):
            return float('inf')

        # Calculate the sum of squared errors if within bounds
        y_pred = growth_curve(A, age, k, B0)
        return np.sum((y - y_pred) ** 2)

    # Call minimize using Nelder-Mead
    result = minimize(loss, initial_params, method=method)
    return result.x

# Cross-validation with option for optimizer choice
def cross_validate(df, optimizer="curve_fit", initial_params=None):
    y = df['agbd'].values
    A = df['nearest_mature_biomass'].values
    age = df['age'].values
    k = df['k_predicted'].values

    if initial_params is None:
        initial_params = [100]
    
    bounds = ([0], 
              [150])

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    r_squared_values = []

    for train_idx, test_idx in kf.split(y):
        k_train, k_test = k[train_idx], k[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        A_train, A_test = A[train_idx], A[test_idx]
        age_train, age_test = age[train_idx], age[test_idx]

        params = fit_minimize(A_train, age_train, k_train, y_train, initial_params, bounds, optimizer)

        # Predict on the test set
        B0, k0 = params
        y_pred = growth_curve(A_test, age_test, k_test, B0)
        print(k0)
        # Calculate R-squared
        ss_res = np.sum((y_test - y_pred) ** 2)
        ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        r_squared_values.append(r_squared)

    return np.mean(r_squared_values)

mean_r_squared_nelder_mead = cross_validate(df, optimizer="Nelder-Mead")
print(f"Mean cross-validated R-squared using Nelder-Mead: {mean_r_squared_nelder_mead}")

# print(np.cumsum(pca.explained_variance_ratio_))
# X_pca_sample, _ = train_test_split(X_pca, test_size=0.8, random_state=42)
# # Step 3: Determine optimal k using Elbow Method and Silhouette Analysis
# inertia_values = []
# silhouette_scores = []
# k_values = range(2, 11)

# for k in k_values:
#     kmeans = KMeans(n_clusters=k, init='k-means++', n_init=20, random_state=42)
#     kmeans.fit(X_pca)
#     inertia_values.append(kmeans.inertia_)
#     silhouette = silhouette_score(X_pca_sample, kmeans.predict(X_pca_sample))
#     silhouette_scores.append(silhouette)
#     print(kmeans.inertia_)
#     print(silhouette)

