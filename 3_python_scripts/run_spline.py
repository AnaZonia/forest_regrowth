import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.interpolate import BSpline
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler

# Function to create spline-transformed features for each predictor
def create_spline_features(X, knots, degree=1):
    n_features = X.shape[1]
    spline_features = []
    for i in range(n_features):
        # Generate a BSpline basis for each feature in X
        coefs = np.zeros(len(knots) - degree - 1)  # Spline coefficient placeholders
        for j in range(len(coefs)):
            coefs[j] = 1.0  # Generate individual basis functions
            spl = BSpline(knots, coefs, k=degree)
            spline_features.append(spl(X[:, i]))
            coefs[j] = 0.0  # Reset for next basis
    return np.column_stack(spline_features)

# Growth curve function
def growth_curve(X, A, age, B0, spline_coefs, theta, knots, degree=2):
    spline_features = create_spline_features(X, knots, degree)
    k = np.dot(spline_features, spline_coefs) * age
    k = np.clip(k, 1e-10, 7)
    result = B0 + (A - B0) * (1 - np.exp(-k)) ** theta
    return result

# Optimization wrapper
def fit_minimize(X, A, age, y, initial_params, bounds, method, knots, degree=2):
    n_spline_coefs = len(create_spline_features(X[:1], knots, degree)[0])

    def loss(params):
        B0, *spline_coefs, theta = params

        # Check if parameters are out of bounds
        if not (bounds[0][0] <= B0 <= bounds[1][0]):
            return float('inf')
        # for i, coef in enumerate(spline_coefs):
        #     if not (bounds[0][i + 1] <= coef <= bounds[1][i + 1]):
        #         return float('inf')
        if not (bounds[0][-1] <= theta <= bounds[1][-1]):
            return float('inf')

        # Calculate the sum of squared errors if within bounds
        y_pred = growth_curve(X, A, age, B0, spline_coefs, theta, knots, degree)
        result = np.sum((y - y_pred) ** 2)
        print(result)
        return np.round(result, decimals=1)

    # Ensure initial_params and bounds match the number of spline coefficients
    if len(initial_params) != n_spline_coefs + 2:  # +2 for B0 and theta
        initial_params = [initial_params[0]] + [0] * n_spline_coefs + [initial_params[-1]]
    
    # if len(bounds[0]) != n_spline_coefs + 2:
    #     bounds = ([bounds[0][0]] + [-5] * n_spline_coefs + [bounds[0][-1]],
    #               [bounds[1][0]] + [5] * n_spline_coefs + [bounds[1][-1]])

    # Call minimize using the specified method
    result = minimize(loss, initial_params, method=method, tol=1e-1)
    return result.x

# Cross-validation with option for optimizer choice
def cross_validate(df, predictors, optimizer="curve_fit", initial_params=None):
    X = df[predictors].values
    y = df['agbd'].values
    A = df['nearest_mature_biomass'].values
    age = df['age'].values

    knots = np.linspace(0, 1, 4)
    degree = 1
    n_spline_coefs = len(create_spline_features(X[:1], knots, degree)[0])

    if initial_params is None:
        initial_params = [100] + [0] * n_spline_coefs + [1.0]
    
    bounds = ([0] + [-5] * n_spline_coefs + [0], 
              [150] + [5] * n_spline_coefs + [10])

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
        params = fit_minimize(X_train, A_train, age_train, y_train, initial_params, bounds, optimizer, knots, degree)

        # Predict on the test set
        B0, *spline_coefs, theta = params
        y_pred = growth_curve(X_test, A_test, age_test, B0, spline_coefs, theta, knots, degree)

        # Calculate R-squared
        ss_res = np.sum((y_test - y_pred) ** 2)
        ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        print(r_squared)
        r_squared_values.append(r_squared)

    return np.mean(r_squared_values)

# Load the data
df = pd.read_csv("./0_data/non_aggregated_all.csv")
df = pd.get_dummies(df, columns=['topography', 'ecoreg', 'indig', 'protec', 'last_LU'], drop_first=True)
df = df.sample(20000, random_state=42)
predictors = ["nitro", "num_fires_before_regrowth", "sur_cover"]

# predictors = [col for col in df.columns if col not in ['age', 'agbd', 'nearest_mature_biomass']]

mean_r_squared_LBFGSB = cross_validate(df, predictors, optimizer="L-BFGS-B")
print(f"Mean cross-validated R-squared using L-BFGS-B: {mean_r_squared_LBFGSB}")
