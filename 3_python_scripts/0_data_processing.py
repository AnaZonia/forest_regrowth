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

    return X, y