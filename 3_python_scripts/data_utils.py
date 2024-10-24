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

# def load_and_preprocess_data(
#         filepath: str, 
#         pars = None, 
#         biome = "both",
#         keep_all_data: bool = False,
#         use_stratified_sample: bool = False,
#         first_stage_sample_size: int = 500,
#         final_sample_size: int = 10000,
#         unseen_portion: float = 0.1
#     ) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, Optional[DataSet]]:
#     """
#     Load and preprocess data from a CSV file.

#     Args:
#         filepath (str): Path to the CSV file.
#         pars (List[str]): List of parameter names to keep.
#         keep_all_data (bool): Flag to keep all data or subset by 'biome'. Defaults to False.
#         use_stratified_sample (bool): Flag to use stratified sampling. Defaults to False.
#         first_stage_sample_size (int): Number of samples per stratification group in the first stage. Defaults to 500.
#         final_sample_size (int): Total sample size to use after stratified sampling. Defaults to 10000.

#     Returns:
#         Tuple[pd.DataFrame, np.ndarray, np.ndarray, Optional[DataSet]]: 
#         Features (X), target (y), asymptote (A), and unseen dataset (if applicable).
#     """
#     df = pd.read_csv(filepath)        
    
#     if biome != "both":
#         df = df[df['biome'] == biome]
#     if biome == 4:
#         df = df.drop(columns = ['mean_aet', 'cwd']) # multicollinearity
#     if biome == "both":
#         df = df.drop(columns = ['mean_aet']) # multicollinearity

#     # Convert 'topography' and 'ecoreg' to categorical if present
#     for col in ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']:
#         if col in df.columns:
#             df[col] = df[col].astype('category')

#     if pars is None:
#         pars = df.columns.drop(["biome", "agbd"]).tolist()

#     if keep_all_data:
#         X = df[pars + ['biome', 'agbd']]
#         unseen_data = None
    
#     if use_stratified_sample:

#         df, unseen_df, pars = stratified_sample_df(df, pars,
#                                                     first_stage_sample_size, final_sample_size, unseen_portion)
#     else:
#         # Split data: 10k rows for df and 1k for unseen_df
#         df, unseen_df = train_test_split(df, test_size = unseen_portion, random_state = 42)

#         df = df.head(final_sample_size)  # Take exactly 10k rows (if needed)
#         unseen_df = unseen_df.head(int(final_sample_size * unseen_portion))  # Take exactly 1k rows

#     X = df[pars]

#     unseen_data = DataSet(
#         X = unseen_df[pars],
#         y = unseen_df['agbd'].values,
#         A = unseen_df['nearest_mature_biomass'].values
#     )

#     y = df['agbd'].values
#     A = df['nearest_mature_biomass'].values  # asymptote

#     return X, y, A, unseen_data


# def stratified_sample_df(df, pars, first_stage_sample_size, final_sample_size, unseen_portion):
#     """
#     Perform stratified sampling on a DataFrame.

#     Args:
#         df (pd.DataFrame): The input DataFrame.
#         pars (List[str]): List of parameter names to keep.
#         test_size (float): Proportion of the data to keep as "unseen".
#         first_stage_sample_size (int): Number of samples per stratification group in the first stage.
#         final_sample_size (int): Total sample size to use after stratified sampling.

#     Returns:
#         Tuple[pd.DataFrame, pd.DataFrame, List[str]]: 
#         Sampled DataFrame, unseen DataFrame, and updated list of parameters.
#     """
#     # Count non-zero values in each column
#     non_zero_counts = (df[pars] != 0).sum()
    
#     # Filter columns based on the count of non-zero values
#     valid_columns = non_zero_counts[non_zero_counts >= first_stage_sample_size].index.tolist()
    
#     # Update pars list
#     pars = [col for col in pars if col in valid_columns]
#     df = df[valid_columns + ['agbd']]

#     # Create a composite stratification variable
#     df['strat_var'] = (df['num_fires_before_regrowth'] > 0).astype(int)
#     # Add contribution from each lulc_sum column
#     lulc_sum_columns = [col for col in df.columns if 'lulc_sum' in col]
#     # Get the top N most common lulc_sum columns
#     top_lulc_columns = df[lulc_sum_columns].sum().nlargest(2).index.tolist()

#     # Add contribution from each of the top lulc_sum columns
#     for i, col in enumerate(top_lulc_columns):
#         df['strat_var'] += (df[col] > 0).astype(int) * (10 ** (i + 1))

#     # Keep 10% of the data as "unseen" for final testing of model performance
#     df, unseen_df = train_test_split(
#         df, test_size = unseen_portion, stratify = df['strat_var'], random_state = 42
#     )

#     unseen_df = unseen_df.head(int(final_sample_size * unseen_portion))

#     # Two-stage sampling
#     # First, ensure all categories are represented
#     sample = df.groupby('strat_var').apply(
#         lambda x: x.sample(min(len(x), first_stage_sample_size), random_state = 42)
#     )

#     # Then, fill the rest with stratified sampling
#     remaining_sample_size = final_sample_size - len(sample)
#     if remaining_sample_size > 0:
#         remaining_df = df[~df.index.isin(sample.index)]
#         stratified_split = StratifiedShuffleSplit(
#             n_splits = 1,
#             test_size = remaining_sample_size,
#             random_state = 42
#         )
#         _, sample_idx = next(stratified_split.split(remaining_df, remaining_df['strat_var']))
#         df = pd.concat([sample, remaining_df.iloc[sample_idx]])
#     return df, unseen_df, pars

def load_and_preprocess_data(
        filepath: str, 
        pars = None, 
        biome = "both",
        ML = False,
        keep_biome: bool = False,
        use_stratified_sample: bool = False,
        first_stage_sample_size: int = 500,
        final_sample_size: int = 10000,
        unseen_portion: float = 0.1
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

    if use_stratified_sample:
        df = df.drop(columns = ['last_LU']) # multicollinearity

    if biome != "both":
        df = df[df['biome'] == biome]
    if biome == 4:
        df = df.drop(columns = ['mean_aet', 'cwd']) # multicollinearity
    if biome == "both":
        df = df.drop(columns = ['mean_aet']) # multicollinearity

    # Select parameters
    if pars is None:
        exclude_cols = ["biome", "agbd"] + ([] if ML else ["nearest_mature_biomass"])
        pars = df.columns.drop(exclude_cols).tolist()

    if keep_biome:
        X = df[pars + ['biome']]
        unseen_data = None

    if use_stratified_sample:
        df, unseen_df, pars = stratified_sample_df(df, pars,
                                                    first_stage_sample_size,
                                                    final_sample_size,
                                                    unseen_portion, ML)
    else:
        # Convert categorical if present
        for col in ['topography', 'ecoreg', 'indig', 'protec']:
            if col in df.columns:
                counts = df[col].value_counts()
                # Identify categories with fewer than 500 occurrences
                rare_categories = counts[counts < first_stage_sample_size].index
                # Remove rows with those categories
                df = df[~df[col].isin(rare_categories)]
                # Convert to categorical
                df[col] = df[col].astype('category')

        # Split data: 10k rows for df and 1k for unseen_df
        df, unseen_df = train_test_split(df, test_size = unseen_portion, random_state = 42)

        df = df.head(final_sample_size)  # Take exactly 10k rows (if needed)
        unseen_df = unseen_df.head(int(final_sample_size * unseen_portion))  # Take exactly 1k rows

    # Prepare features, target, and asymptote
    X = df[pars]
    y = df['agbd'].values
    A = df['nearest_mature_biomass'].values

    # Prepare unseen dataset
    unseen_data = DataSet(
        X=unseen_df[pars],
        y=unseen_df['agbd'].values,
        A=unseen_df['nearest_mature_biomass'].values
    )

    return X, y, A, unseen_data



# Stratified sampling function
def stratified_sample_df(df, pars, first_stage_sample_size, final_sample_size, unseen_portion, ML):
    """
    Perform stratified sampling on a DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.
        pars (List[str]): List of parameter names to keep.
        first_stage_sample_size (int): Number of samples per stratification group in the first stage.
        final_sample_size (int): Total sample size to use after stratified sampling.
        unseen_portion (float): Proportion of data to keep as "unseen" for testing.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, List[str]]: 
        Sampled DataFrame, unseen DataFrame, and updated list of parameters.
    """
    # Add unique identifier
    df['original_index'] = df.index
    # Store original relationships for verification
    original_df = df.copy()

    # Remove categories with fewer than 100 occurrences
    categorical_colnames = ['topography', 'ecoreg', 'indig', 'protec']
    for col in categorical_colnames:
        if col in df.columns:
            counts = df[col].value_counts()
            rare_categories = counts[counts < 100].index
            df = df[~df[col].isin(rare_categories)]
    
    # Remove columns with fewer than 500 non-zero values
    non_zero_counts = (df[pars] != 0).sum()
    valid_columns = non_zero_counts[non_zero_counts >= first_stage_sample_size].index.tolist()

    # Update pars list
    pars = [col for col in pars if col in valid_columns]
    keep_cols = ['original_index', 'agbd'] + ([] if ML else ["nearest_mature_biomass"])
    df = df[keep_cols + valid_columns].copy()

    # Create a composite stratification variable
    df.loc[:, 'strat_var'] = (df['num_fires_before_regrowth'] > 0).astype(int)
    # Add contribution from each lulc_sum column and each category
    lulc_sum_columns = [col for col in df.columns if 'lulc_sum' in col]
    selected_columns = lulc_sum_columns #+ [col for col in categorical_colnames if col in df.columns]

    # Add contribution from each of the selected_columns
    for i, col in enumerate(selected_columns):
        df.loc[:, 'strat_var'] += (df[col] > 0).astype(int) * (10 ** (i + 1))

    # Remove rows with infrequent strat_var values
    value_counts = df['strat_var'].value_counts()
    values_to_remove = value_counts[value_counts < 50].index
    df = df[~df['strat_var'].isin(values_to_remove)]

    # Convert categorical columns to category type
    for col in categorical_colnames:
        if col in df.columns:
            df[col] = df[col].astype('category')

    # Keep 10% of the data as "unseen" for final testing of model performance
    df, unseen_df = train_test_split(
        df, test_size = unseen_portion, stratify = df['strat_var'], random_state = 42
    )
    unseen_df = unseen_df.head(int(final_sample_size * unseen_portion))

    # Two-stage sampling
    # First, ensure all categories are represented
    sample = df.groupby('strat_var').apply(
        lambda x: x.sample(min(len(x), first_stage_sample_size), random_state = 42)
    )

    # Then, fill the rest with stratified sampling
    remaining_sample_size = final_sample_size - len(sample)
    if remaining_sample_size > 0:
        # Exclude rows that have already been sampled in the first stage
        remaining_df = df[~df['original_index'].isin(sample['original_index'])].reset_index(drop=True)

        stratified_split = StratifiedShuffleSplit(
            n_splits=1,
            test_size = remaining_sample_size,  # Oversample by 50%
            random_state=42
        )

        _, indices = next(stratified_split.split(remaining_df, remaining_df['strat_var']))
        
        # Concatenate the initial sample and the newly sampled data, avoiding duplicates
        df = pd.concat([sample, remaining_df.loc[indices]])

    verify_data_integrity(original_df, df, unseen_df, pars)
    df = df.drop(columns = ['strat_var', 'original_index'])
    unseen_df = unseen_df.drop(columns = ['strat_var', 'original_index'])

    return df, unseen_df, pars

# Data integrity verification function
def verify_data_integrity(original_df, sampled_df, unseen_df, pars):
    """
    Verify the integrity of the sampled and unseen data.

    Args:
        original_df (pd.DataFrame): Original DataFrame.
        sampled_df (pd.DataFrame): Sampled DataFrame.
        unseen_df (pd.DataFrame): Unseen DataFrame.
        pars (List[str]): List of parameter names.

    Raises:
        AssertionError: If any integrity check fails.
    """
    # Check for duplicate indices
    assert len(sampled_df['original_index']) == len(sampled_df['original_index'].unique()), "Duplicate indices found in sampled data"
    assert len(unseen_df['original_index']) == len(unseen_df['original_index'].unique()), "Duplicate indices found in unseen data"

    # Check if all sampled rows exist in original data
    all_sampled_indices = set(sampled_df['original_index']).union(set(unseen_df['original_index']))
    assert all_sampled_indices.issubset(set(original_df['original_index'])), "Some sampled rows don't exist in original data"

    # Verify relationships between columns
    original_relationships = original_df[['original_index', 'agbd'] + pars].head(10).to_dict('records')
    for record in original_relationships:
        orig_index = record['original_index']
        if orig_index in sampled_df['original_index'].values:
            sampled_record = sampled_df[sampled_df['original_index'] == orig_index].iloc[0]
        elif orig_index in unseen_df['original_index'].values:
            sampled_record = unseen_df[unseen_df['original_index'] == orig_index].iloc[0]
        else:
            continue

        for col in pars + ['agbd']:
            if col in sampled_record:
                assert np.isclose(record[col], sampled_record[col]), f"Mismatch in {col} for index {orig_index}"

    print("All integrity checks passed!")


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