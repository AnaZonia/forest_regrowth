# scripts/data_processing.py
import pandas as pd

def load_and_preprocess_data(filepath):
    
    df = pd.read_csv(filepath)
    # Assuming df is your DataFrame
    pars = [
        "age", "nearest_mature", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
        "lulc_sum_40", "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
        "cwd"
    ]

    # Select only the specified columns
    X = df[pars]    
    y = df['agbd']

    # Create the new DataFrame with one row
    pars = pd.DataFrame(columns=pars)
    # Fill the row with values
    pars.loc[0] = 0  # Set all values to 1 initially
    pars.loc[0, 'theta'] = 1
    pars.loc[0, 'B0'] = df['agbd'].mean()

    return X, y, pars