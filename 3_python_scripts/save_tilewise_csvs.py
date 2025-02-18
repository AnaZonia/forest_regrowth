import glob
import os

# Define the folder path
folder_path = './0_data/tilewise_exports'

# Get a list of all CSV files in the folder
all_csv_files = glob.glob(os.path.join(folder_path, '*.csv'))

# Separate lists for "15yr" and "10yr" files
files_5 = [f for f in all_csv_files if '_5yr' in os.path.basename(f)]
files_10 = [f for f in all_csv_files if '_10yr' in os.path.basename(f)]

# Function to read files and exclude all-NA DataFrames
def read_non_na_csv(files):
    dfs = [pd.read_csv(f) for f in files]
    return [df for df in dfs if not df.empty and not df.isna().all(axis=None)]

# Read and concatenate (rbind) non-NA "15yr" files
files_5_dfs = read_non_na_csv(files_5)
if files_5_dfs:
    files_5_df = pd.concat(files_5_dfs, ignore_index=True)
else:
    files_5_df = pd.DataFrame()  # Empty DataFrame if no non-NA files

# Read and concatenate (rbind) non-NA "10yr" files
files_10_dfs = read_non_na_csv(files_10)
if files_10_dfs:
    files_10_df = pd.concat(files_10_dfs, ignore_index=True)
else:
    files_10_df = pd.DataFrame()  # Empty DataFrame if no non-NA files

# Verify the results (optional)
print("files_5_df DataFrame shape:", files_5_df.shape)
print("files_10_df DataFrame shape:", files_10_df.shape)

files_5_df.to_csv('./0_data/aggregated_5yr_tilewise.csv', index=False)
print("Saved aggregated 5yr data to 0_data/aggregated_5yr_tilewise.csv")
files_10_df.to_csv('./0_data/aggregated_10yr_tilewise.csv', index=False)
print("Saved aggregated 10yr data to 0_data/aggregated_10yr_tilewise.csv")
    
