import pandas as pd
import numpy as np
from typing import List, Tuple, NamedTuple, Optional
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit

class DataSet(NamedTuple):
    X: pd.DataFrame
    y: np.ndarray
    A: np.ndarray

def load_and_preprocess_data(
        filepath: str, 
        pars: List[str], 
        keep_all_data: bool = False,
        use_stratified_sample: bool = False,
        first_stage_sample_size: int = 500,
        final_sample_size: int = 10000
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

    if keep_all_data:
        X = df[pars + ['biome', 'nearest_mature']]
        unseen_data = None
    else:
        if use_stratified_sample:
                # Keep 10% of the data as "unseen" for final testing of model performance
            test_size = len(df) * 0.1

            df, unseen_df, pars = stratified_sample_df(df, pars, test_size, \
                                                       first_stage_sample_size, final_sample_size)
        else:
            # df = df[df['biome'] == 1]

            # Split data: 10k rows for df and 1k for unseen_df
            df, unseen_df = train_test_split(df, test_size = 0.1, random_state = 42)

            df = df.head(10000)  # Take exactly 10k rows (if needed)
            unseen_df = unseen_df.head(1000)  # Take exactly 1k rows

        X = df[pars]

        unseen_data = DataSet(
            X = unseen_df[pars],
            y = unseen_df['agbd'].values,
            A = unseen_df['nearest_mature'].values
        )

    y = df['agbd'].values
    A = df['nearest_mature'].values  # asymptote

    return X, y, A, unseen_data


def stratified_sample_df(df, pars, test_size,
        first_stage_sample_size, 
        final_sample_size):

    # Count non-zero values in each column
    non_zero_counts = (df[pars] != 0).sum()
    
    # Filter columns based on the count of non-zero values
    valid_columns = non_zero_counts[non_zero_counts >= first_stage_sample_size].index.tolist()
    
    # Update pars list
    pars = [col for col in pars if col in valid_columns]

    # Create a composite stratification variable
    df['strat_var'] = (
        (df['lulc_sum_39'] > 0).astype(int) * 100 +
        (df['num_fires_before_regrowth'] > 0).astype(int) * 10 +
        (df['lulc_sum_41'] > 0).astype(int)
    )

    # Keep 10% of the data as "unseen" for final testing of model performance
    df, unseen_df = train_test_split(
        df, test_size = test_size, stratify = df['strat_var'], random_state = 42
    )

    # Two-stage sampling
    # First, ensure all categories are represented
    sample = df.groupby('strat_var').apply(
        lambda x: x.sample(min(len(x), first_stage_sample_size), random_state = 42)
    )

    # Then, fill the rest with stratified sampling
    remaining_sample_size = final_sample_size - len(sample)
    if remaining_sample_size > 0:
        remaining_df = df[~df.index.isin(sample.index)]
        stratified_split = StratifiedShuffleSplit(
            n_splits = 1,
            test_size = remaining_sample_size,
            random_state = 42
        )
        _, sample_idx = next(stratified_split.split(remaining_df, remaining_df['strat_var']))
        df = pd.concat([sample, remaining_df.iloc[sample_idx]])
    return df, unseen_df, pars


def make_initial_parameters(pars, y, lag = False):
    if lag:
        initial_params = np.zeros(len(pars) + 3)
        initial_params[0] = 0  # m_base
        initial_params[1] = 1  # sd_base
        initial_params[2] = 1  # sd
        initial_params[3] = 1  # theta
        initial_params[4] = 0.0001  # par1
    else:
        initial_params = np.zeros(len(pars) + 2)
        initial_params[0] = y.mean()  # B0
        initial_params[1] = 1  # theta

    return initial_params

def format_best_params(best_params, pars, func_form):

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
            *[best_params[f"coeff_{i}"] for i in range(len(pars))]
        ])
    else:
        # Fallback: treat all as coefficients if no special params are found
        params = np.array([
            *[best_params[f"coeff_{i}"] for i in range(len(pars))]
        ])

    return params