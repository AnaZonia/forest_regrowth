
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, cross_validate, GridSearchCV, learning_curve
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit


def load_and_preprocess_data(
        filepath: str, 
        pars = None, 
        keep_all_data: bool = False,
        final_sample_size: int = 10000
    ) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Load and preprocess data from a CSV file.

    Args:
        filepath (str): Path to the CSV file.
        pars (List[str]): List of parameter names to keep.
        keep_all_data (bool): Flag to keep all data or subset by 'biome'. Defaults to False.
        final_sample_size (int): Total sample size to use after stratified sampling. Defaults to 10000.

    Returns:
        Tuple[pd.DataFrame, np.ndarray]: Features (X) and target (y).
    """
    df = pd.read_csv(filepath)        

    df = df.loc[:, ~df.columns.str.contains(r"\d{4}")]
    
    # Identify columns with fewer than 200 non-zero values
    low_value_columns = [col for col in df.columns if (df[col] != 0).sum() < 100]
    df = df.drop(columns = low_value_columns)

    # Convert 'topography' and 'ecoreg' to categorical if present
    for col in ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']:
        if col in df.columns:
            df[col] = df[col].astype('category')
            # Filter out categories with fewer than 50 rows
            counts = df[col].value_counts()
            valid_categories = counts[counts >= 50].index
            df = df[df[col].isin(valid_categories)]

    if pars is None:
        pars = df.columns.tolist()
    
    if not keep_all_data:
        pars = [col for col in pars if col not in ["biome", "biomass", "latitude", "longitude", "silva_age", "heinrich_biomass_2020", 'amazon_', 'last_LU', 'fallow', "regions", "mean_biomass_quarter", "Unnamed: 0", "X"]]

    X = df[pars]
    y = df['biomass'].values
    # lat_lon = df[['latitude', 'longitude']]

    return X, y#, lat_lon



def print_feature_importance(perm_importance, feature_names):
    """
    Print the feature importance based on permutation importance.

    Parameters:
    - perm_importance: PermutationImportance object containing importances.
    - feature_names: List of feature names.
    """
    feature_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': perm_importance.importances_mean,
        'std': perm_importance.importances_std
    })
    feature_importance = feature_importance.sort_values('importance', ascending = False)

    print("\nFeature Importance:")
    print(feature_importance.to_string(index = False))
    plot_feature_importance(feature_importance)




def plot_feature_importance(feature_importance):
    """
    Plot feature importance with error bars.

    Parameters:
    - feature_importance: DataFrame containing features, importances, and std deviations.
    """
    plt.figure(figsize=(10, 6))
    plt.barh(feature_importance['feature'], feature_importance['importance'], 
             xerr=feature_importance['std'], color='skyblue', edgecolor='black')
    
    plt.xlabel('Importance')
    plt.title('Feature Importance with Standard Deviation')
    plt.gca().invert_yaxis()  # Invert y axis to have the most important feature on top
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()



def cross_validation(X, y, model, param_grid = None):
    """
    Perform cross-validation for regression models.

    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        model: The regression model to train.
        param_grid (dict): Parameter grid for hyperparameter tuning.

    Returns:
        tuple: Mean R2, standard deviation of R2, permutation importance, predictions array.
    """
    model = Pipeline([
        ('scaler', MinMaxScaler()),
        ('regressor', model)
    ])

    kf = KFold(n_splits = 5, shuffle = True, random_state = 42)
    splits = list(kf.split(X))
    
    if param_grid:
        param_grid = {'regressor__' + k: v for k, v in param_grid.items()}
        grid_search = GridSearchCV(model, param_grid, cv = kf, 
                            scoring = 'neg_mean_squared_error', refit = True)
        grid_search.fit(X, y)

        best_params = grid_search.best_params_
        model = model.set_params(**best_params)
    else:
        model = model.fit(X, y)


    # After performing grid search, perform cross-validation
    cv_results = cross_validate(model, X, y, cv = kf, 
                        scoring = ['neg_mean_squared_error'],
                        return_estimator = True)
    mse_scores = -cv_results['test_neg_mean_squared_error']

    # Perform permutation importance
    best_index = np.argmin(-cv_results['test_neg_mean_squared_error'])
    _, test_indices = splits[best_index]
    X_test = X.iloc[test_indices]
    y_test = y[test_indices]
    best_model = cv_results['estimator'][best_index]
    perm_importance = permutation_importance(best_model, X_test, y_test, n_repeats = 10, random_state = 42)

    # Calculate R2 from MSE
    r2_scores = 1 - (mse_scores / np.var(y))
    r2_mean = np.mean(r2_scores)
    r2_sd = np.std(r2_scores)

    predictions = np.empty(len(X))
    predictions[:] = np.nan
    for i, estimator in enumerate(cv_results['estimator']):
        _, test_indices = list(kf.split(X))[i]
        fold_predictions = estimator.predict(X.iloc[test_indices])
        predictions[test_indices] = fold_predictions

    # Calculate R-squared between all biomass values (y) and all predictions across folds
    r2_final = r2_score(y, predictions)

    return r2_mean, r2_sd, r2_final, perm_importance, predictions


def run_model(datasource, models, param_grids, filepath, pars = None):

    # filepath = f"./0_data/{datasource}.csv"

    # Load and preprocess data for the given biome
    X, y = load_and_preprocess_data(filepath, pars, \
                        final_sample_size = 8000)

    rows = pd.DataFrame()
    # Perform regression analysis for each model
    for name, model in models.items():
        r2_mean, r2_sd, r2_final, perm_importance, _ = regression_cv(
            X, y, model, param_grids[name]
        )
        
        print_feature_importance(perm_importance, X.columns)

        row = {
            'datasource': datasource,
            'model': name,
            'r2_mean': r2_mean,
            'r2_sd': r2_sd,
            'r2_final': r2_final
        }

        print(row)
        # Append result to the DataFrame
        rows = pd.concat([rows, pd.DataFrame([row])], ignore_index=True)

    return rows




if __name__ == "__main__":

    models = {
        "Linear Regression": LinearRegression(),
        "XGBoost": XGBRegressor(random_state = 42),
        "Random Forest": RandomForestRegressor(random_state = 42)
    }

    param_grids = {
        "Linear Regression": None,
        "XGBoost": {
            'learning_rate': [0.01, 0.1],
            'n_estimators': [100, 500],
            'max_depth': [3, 7],
            'min_child_weight': [1, 5]
        },
        "Random Forest": {
            'n_estimators': [100, 300],
            'max_depth': [None, 20],
            'min_samples_split': [2, 10],
            'min_samples_leaf': [1, 4]
        }
    }


    # # Define keywords for land use and landscape parameters, and include "age"
    # land_use_keywords = ['lulc', 'LU', 'fallow', 'num_fires']
    # landscape_keywords = ['sur_cover', 'distance', 'nearest_mature_biomass']
    # essential_keywords = land_use_keywords + landscape_keywords + ['age']

    results_df = pd.DataFrame()

    datasource = "~/Documents/data/mapbiomas_heinrich_lulc_regions_renamed.csv"

    # pars = ["age", "cwd", "phh2o", "mean_pr", "mean_srad", "mean_temp", "mean_aet", "mean_vpd", "mean_soil", "nearest_mature_biomass"]
    # pars = ["age", "cwd", "phh2o", "mean_pr", "mean_srad", "mean_temp", "mean_aet", "mean_vpd", "mean_soil", "nearest_mature_biomass",
    # "protec", "indig", "nitro", "sand", "ecoreg", "topography"]
    # pars = ["age", "cwd", "phh2o", "mean_pr", "mean_srad", "mean_temp", "mean_aet", "mean_vpd", "mean_soil", "nearest_mature_biomass",
    # "protec", "indig", "nitro", "sand", "ecoreg", "topography", "distance", "sur_cover"]
    # pars = ["age", "cwd", "phh2o", "mean_pr", "mean_srad", "mean_temp", "mean_aet", "mean_vpd", "mean_soil", "nearest_mature_biomass",
    # "protec", "indig", "nitro", "sand", "ecoreg", "topography", "distance", "sur_cover", "num_fires"]
    pars = ["age", "cwd", "phh2o", "mean_pr", "mean_srad", "mean_temp", "mean_aet", "mean_vpd", "mean_soil", "nearest_mature_biomass",
    "protec", "indig", "nitro", "sand", "ecoreg", "topography", "distance", "sur_cover", "num_fires",
    "lulc_sum_10", "lulc_sum_30", "last_LU", "fallow"]

 

    results = run_model(
        datasource,
        {"XGBoost": models["XGBoost"]},
        {"XGBoost": param_grids["XGBoost"]},
        datasource,
        pars
    )

    # Append results and update progress
    results_df = pd.concat([results_df, results], ignore_index=True)

    # Save CSV after each iteration
    results_df.to_csv("./0_results/feb2nd_results.csv", index=False)
    print("Data saved to CSV.")

