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
        test_size: int = 10000, 
        first_stage_sample_size: int = 500, 
        final_sample_size: int = 10000
    ) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, Optional[DataSet]]:
    """
    Load and preprocess data from a CSV file.

    Args:
        filepath (str): Path to the CSV file.
        pars (List[str]): List of parameter names to keep.
        keep_all_data (bool): Flag to keep all data or subset by 'biome'. Defaults to False.
        test_size (int): Number of samples to use as unseen data. Defaults to 10000.
        first_stage_sample_size (int): Number of samples per stratification group in the first stage. Defaults to 100.
        final_sample_size (int): Total sample size to use after stratified sampling. Defaults to 10000.

    Returns:
        Tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, Optional[DataSet]]: 
        Features (X), target (y), asymptote (A), initial parameters, and unseen dataset (if applicable).
    """
    df = pd.read_csv(filepath)        

    if keep_all_data:
        X = df[pars + ['biome', 'nearest_mature']]
        unseen_data = None
    else:
        df = df[df['biome'] == 1]
        # df = df[df['num_fires_after_regrowth'] == 0]

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
            df, test_size=test_size, stratify=df['strat_var'], random_state=42
        )

        # Two-stage sampling
        # First, ensure all categories are represented
        sample = df.groupby('strat_var').apply(
            lambda x: x.sample(min(len(x), first_stage_sample_size), random_state=42)
        )

        # Then, fill the rest with stratified sampling
        remaining_sample_size = final_sample_size - len(sample)
        if remaining_sample_size > 0:
            remaining_df = df[~df.index.isin(sample.index)]
            stratified_split = StratifiedShuffleSplit(
                n_splits=1,
                test_size=remaining_sample_size,
                random_state=42
            )
            _, sample_idx = next(stratified_split.split(remaining_df, remaining_df['strat_var']))
            df = pd.concat([sample, remaining_df.iloc[sample_idx]])

        X = df[pars]

        unseen_data = DataSet(
            X=unseen_df[pars],
            y=unseen_df['agbd'].values,
            A=unseen_df['nearest_mature'].values
        )
    

    y = df['agbd'].values
    A = df['nearest_mature'].values  # asymptote
    print((X != 0).sum())

    initial_params = np.zeros(len(pars) + 2)
    initial_params[0] = y.mean()
    initial_params[1] = 1

    return X, y, initial_params, A, unseen_data


def load_and_preprocess_data_simplified(
        filepath: str, 
        pars: List[str], 
        test_size: int = 10000 
    ):

    df = pd.read_csv(filepath)
    df = df.rename(columns={'age_eu': 'age'})

    # Keep 10% of the data as "unseen" for final testing of model performance
    df, unseen_df = train_test_split(
        df, test_size=test_size, random_state=42
    )

    X = df[pars]

    unseen_data = DataSet(
        X=unseen_df[pars],
        y=unseen_df['agbd'].values,
        A=unseen_df['nearest_mature'].values
    )

    y = df['agbd'].values
    A = df['nearest_mature'].values  # asymptote

    initial_params = np.zeros(len(pars) + 2)
    initial_params[0] = y.mean() # B0
    initial_params[1] = 1 # theta

    return X, y, initial_params, A, unseen_data