"""
Nelder-Mead Optimization and Regression Analysis Script

This script performs regression analysis and Nelder-Mead optimization on a dataset.
It includes functionality for:
1. Standard regression analysis using Linear Regression, XGBoost, and Random Forest
2. Nelder-Mead optimization with optional hyperparameter tuning
3. Visualization of results

Modules:
- data_utils: Contains utility functions for loading, preprocessing, and handling data.
- model_utils: Contains functions for performing regression and Nelder-Mead optimization.
- tuners: Contains functions for hyperparameter tuning using various methods.

Functions:
- print_feature_importance: Prints the feature importance based on permutation importance.
- regression_main: Main function to perform regression analysis using different models.
- nelder_mead_main: Main function to perform Nelder-Mead optimization with optional parameter tuning.

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

import pandas as pd
from typing import Tuple, Optional
from sklearn.model_selection import train_test_split


from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, cross_validate, GridSearchCV, learning_curve
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance


# DataSet class definition
class DataSet:
    def __init__(self, X, y, A):
        self.X = X  # Features
        self.y = y  # Target variable
        self.A = A  # Asymptote

def load_and_preprocess_data(
        filepath: str, 
        pars = None, 
        biome = "both",
        keep_all_data: bool = False,
        final_sample_size: int = 10000,
        unseen_portion: float = 0.1,
        remove_land_use_params: bool = False  # New switch
    ) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, Optional[DataSet]]:
    """
    Load and preprocess data from a CSV file.

    Args:
        filepath (str): Path to the CSV file.
        pars (List[str]): List of parameter names to keep.
        keep_all_data (bool): Flag to keep all data or subset by 'biome'. Defaults to False.
        use_stratified_sample (bool): Flag to use stratified sampling. Defaults to False.
        first_stage_sample_size (int): Number of samples per stratification group in the first stage. Defaults to 500.
        final_sample_size (int): Total sample size to use after stratified sampling. Defaults to 10000.

    Returns:
        Tuple[pd.DataFrame, np.ndarray, np.ndarray, Optional[DataSet]]: 
        Features (X), target (y), asymptote (A), and unseen dataset (if applicable).
    """
    df = pd.read_csv(filepath)        
    
    if biome != "both":
        df = df[df['biome'] == biome]
    # if biome == 4:
    #     df = df.drop(columns = ['mean_aet', 'cwd']) # multicollinearity
    # if biome == "both":
    #     df = df.drop(columns = ['mean_aet']) # multicollinearity

    # Convert 'topography' and 'ecoreg' to categorical if present
    for col in ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']:
        if col in df.columns:
            df[col] = df[col].astype('category')

    if pars is None:
        pars = df.columns.drop(["biome", "agbd"]).tolist()

    if keep_all_data:
        X = df[pars + ['biome', 'agbd']]
        unseen_data = None

    # Remove land use related parameters if the switch is on
    if remove_land_use_params:
        keywords = ['lulc', 'LU', 'fallow', 'mature_forest_years', 'num_fires', 'cover']
        pars = [par for par in pars if not any(keyword in par for keyword in keywords)]

    # Split data: 10k rows for df and 1k for unseen_df
    df, unseen_df = train_test_split(df, test_size = unseen_portion, random_state = 42)

    df = df.sample(final_sample_size, random_state = 42)  # Take exactly 10k rows (if needed)
    unseen_df = unseen_df.sample(n = int(final_sample_size * unseen_portion), random_state = 42)  # Take exactly 1k rows

    X = df[pars]

    unseen_data = DataSet(
        X = unseen_df[pars],
        y = unseen_df['agbd'].values,
        A = unseen_df['nearest_mature_biomass'].values
    )

    y = df['agbd'].values
    A = df['nearest_mature_biomass'].values  # asymptote

    return X, y, A, unseen_data



def plot_learning_curves(X, y, model, name, cv = 5):
    """
    Plot learning curves for a given model pipeline.
    
    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        model: The model to train.
        name (str): Name of the model for the plot title.
        cv (int): Number of cross-validation folds.
    """
    train_sizes, train_scores, test_scores = learning_curve(
        model, X, y, cv = cv, scoring = 'neg_mean_squared_error',
        train_sizes = np.linspace(0.1, 1.0, 10), n_jobs = -1)
    
    # Calculate mean and standard deviation for training set scores
    train_mean = -np.mean(train_scores, axis = 1)
    train_std = np.std(train_scores, axis = 1)

    # Calculate mean and standard deviation for test set scores
    test_mean = -np.mean(test_scores, axis = 1)
    test_std = np.std(test_scores, axis = 1)

    # Plot learning curve
    plt.figure(figsize = (10, 6))
    plt.fill_between(train_sizes, train_mean - train_std, train_mean + train_std, alpha = 0.1, color = "r")
    plt.fill_between(train_sizes, test_mean - test_std, test_mean + test_std, alpha = 0.1, color = "g")
    plt.plot(train_sizes, train_mean, 'o-', color = "r", label = "Training score")
    plt.plot(train_sizes, test_mean, 'o-', color = "g", label = "Cross-validation score")
    
    plt.xlabel("Training examples")
    plt.ylabel("Mean Squared Error")
    plt.title(f"Learning Curves for {name}")
    plt.legend(loc = "best")
    plt.grid(True)
    
    return plt



def regression_cv(X, y, model, unseen_data, name, param_grid = None):
    """
    Perform cross-validation for regression models.

    Args:
        X (pd.DataFrame): Feature matrix.
        y (np.ndarray): Target values.
        model: The regression model to train.
        unseen_data (DataSet): Unseen data for evaluation.
        name (str): Name of the model for the plot title.
        param_grid (dict): Parameter grid for hyperparameter tuning.

    Returns:
        tuple: Mean R2, standard deviation of R2, unseen R2, permutation importance, and learning curve figure.
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
                            scoring = 'neg_mean_squared_error', refit = False)
        grid_search.fit(X, y)

        best_params = grid_search.best_params_
        model = model.set_params(**best_params)

    # After performing grid search, perform cross-validation
    cv_results = cross_validate(model, X, y, cv = kf, 
                        scoring = ['neg_mean_squared_error'],
                        return_estimator = True)
    mse_scores = -cv_results['test_neg_mean_squared_error']
    best_index = np.argmin(-cv_results['test_neg_mean_squared_error'])
    best_model = cv_results['estimator'][best_index]

    # Perform permutation importance
    _, test_indices = splits[best_index]
    X_test = X.iloc[test_indices]
    y_test = y[test_indices]
    perm_importance = permutation_importance(best_model, X_test, y_test, n_repeats = 10, random_state = 42)

    # Calculate R2 from MSE
    r2_scores = 1 - (mse_scores / np.var(y))
    mean_r2 = np.mean(r2_scores)
    std_r2 = np.std(r2_scores)

    # Calculate r2 of the best model in the unseen dataset
    y_pred = best_model.predict(unseen_data.X)
    unseen_r2 = r2_score(unseen_data.y, y_pred)

    fig = plot_learning_curves(X, y, best_model, name)

    return mean_r2, std_r2, unseen_r2, perm_importance, fig, best_model


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


def regression_main():
    """
    Main function to perform regression analysis using different models for multiple biomes and data sources.
    Saves the results into a CSV file.
    """
    # Define biomes and data sources
    biomes = [1, 4, "both"]  # You can extend this as needed
    filepaths = {
        # "eu": "./0_data/eu.csv",
        "mapbiomas": "./0_data/non_aggregated.csv"
    }

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

    results_list = []
    # Loop over each biome and data source
    for biome in biomes:
        for datasource, filepath in filepaths.items():
            # Load and preprocess data for the given biome
            X, y, _, unseen_data = load_and_preprocess_data(filepath, \
                                biome = biome,
                                final_sample_size = 15000,
                                unseen_portion = 0.2,
                                remove_land_use_params = True)

            # Perform regression analysis for each model
            for name, model in models.items():
                mean_r2, std_r2, unseen_r2, perm_importance, _, best_model = regression_cv(
                    X, y, model, unseen_data, name, param_grids[name]
                )
                
                # Print results for this combination
                print(f"\n{name} Results for Biome {biome}, DataSource: {datasource}")
                print(f"Cross-validation R2: {mean_r2:.3f} (±{std_r2:.3f})")
                print(f"Unseen data R2: {unseen_r2:.3f}")
                
                print_feature_importance(perm_importance, X.columns)

                results_list.append({
                    'biome': biome,
                    'datasource': datasource,
                    'model': name,
                    'cv_r2': mean_r2,
                    'unseen_r2': unseen_r2
                })

                # Add predictions to the original dataset
                X['pred'] = best_model.predict(X)
                X['agbd'] = y
                # Save the dataset with predictions
                output_filepath = f"./0_results/predictions_{datasource}_biome_{biome}.csv"
                X.to_csv(output_filepath, index=False)
                print(f"Predictions saved to {output_filepath}")

                                
                # # Display the figure (optional)
                # fig.suptitle(f"{name} - Biome {biome} - {datasource}")
                # plt.show()

    # # After the loop, create the DataFrame from the list of results
    # results_df = pd.DataFrame(results_list)
    # # Save results to CSV
    # results_df.to_csv("./0_results/regression_results_nocat.csv", index = False)
    # print("Results saved to ./0_results/regression_results.csv")




if __name__ == "__main__":
    regression_main()