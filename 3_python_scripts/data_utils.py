"""
data_utils.py

This module contains utility functions for loading, preprocessing, and handling data for regression analysis.
It includes functions for loading data from a CSV file, performing stratified sampling, and generating initial
parameters for optimization algorithms.

Functions:
- load_and_preprocess_data: Load and preprocess data from a CSV file.
- stratified_sample_df: Perform stratified sampling on a DataFrame.
- make_initial_parameters: Generate initial parameters for optimization.
- format_best_params: Format the best parameters for optimization.

Classes:
- DataSet: NamedTuple to store features, target, and asymptote data.
"""

import pandas as pd
import numpy as np
from typing import Tuple, Optional, List
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit


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
        # exclude_cols = ["biome", "agbd"] + ([] if ML else ["nearest_mature_biomass"])
        # pars = df.columns.drop(exclude_cols).tolist()

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




def make_initial_parameters(pars, y, A, func_form):
    """
    Generate initial parameters for optimization.

    Args:
        pars (List[str]): List of parameter names.
        y (np.ndarray): Target values.
        lag (bool): Flag to include lag parameters. Defaults to False.

    Returns:
        np.ndarray: Initial parameters for optimization.
    """
    if func_form == "lag":
        initial_params = np.full(len(pars) + 4, 0.0001)
        initial_params[0] = 0  # m_base
        initial_params[1] = 1  # sd_base
        initial_params[2] = 1  # sd
        initial_params[3] = 1  # theta
        # initial_params[4] = -np.log(1 - (y.mean() / A.mean())) # k0
    else:
        initial_params = np.zeros(len(pars) + 2)
        initial_params[0] = y.mean()  # B0
        initial_params[1] = 1  # theta

    return initial_params

def format_best_params(best_params, pars, func_form):
    """
    Format the best parameters for optimization.

    Args:
        best_params (dict): Dictionary of best parameters.
        pars (List[str]): List of parameter names.
        func_form (str): Functional form ('B0_theta' or 'lag').

    Returns:
        np.ndarray: Formatted parameters for optimization.
    """
    if func_form == "B0_theta":
        params = np.array([
            best_params["B0"],
            best_params["theta"],
            *[best_params[f"coeff_{i}"] for i in range(len(pars))]
        ])
    elif func_form == "lag":
        params = np.array([
            best_params["m_base"],
            best_params["sd_base"],
            best_params["sd"],
            best_params["theta"],
            *[best_params[f"coeff_{i}"] for i in range(len(pars) - 1)]
        ])
    else:
        # Fallback: treat all as coefficients if no special params are found
        params = np.array([
            *[best_params[f"coeff_{i}"] for i in range(len(pars))]
        ])

    return params