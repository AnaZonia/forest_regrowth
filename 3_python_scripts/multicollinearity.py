import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.outliers_influence import variance_inflation_factor

# Import the load_and_preprocess_data function from your existing script
from fit_XGBoost import load_and_preprocess_data

def calculate_correlation_matrix(df):
    """Calculate the correlation matrix for the dataframe."""
    return df.corr()

def plot_correlation_heatmap(corr_matrix):
    """Plot a heatmap of the correlation matrix."""
    plt.figure(figsize = (12, 10))
    sns.heatmap(corr_matrix, annot = True, cmap = 'coolwarm', vmin = -1, vmax = 1, center = 0)
    plt.title('Correlation Heatmap')
    plt.tight_layout()
    plt.show()

def calculate_vif(df):
    """Calculate the Variance Inflation Factor for each feature."""
    scaler = StandardScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df), columns = df.columns)
    
    vif_data = pd.DataFrame()
    vif_data["feature"] = df.columns
    vif_data["VIF"] = [variance_inflation_factor(df_scaled.values, i) for i in range(df_scaled.shape[1])]
    
    return vif_data.sort_values('VIF', ascending = False)

def main():
    # Define the parameters (features) you want to use
    # pars = [
    #     "nearest_mature", "age", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
    #     "lulc_sum_40", "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
    #     "cwd"
    # ]
    # Load and preprocess the data
    # X, _, _, _ = load_and_preprocess_data("./0_data/eu.csv")#, use_stratified_sample = True)
    X, _, _, _ = load_and_preprocess_data("./0_data/non_aggregated.csv", \
                                           biome = 1, keep_all_data = True,
                                           use_stratified_sample = True,
                                first_stage_sample_size = 500, final_sample_size = 15000,
                                unseen_portion = 0.2)

    # Calculate and plot correlation matrix
    corr_matrix = calculate_correlation_matrix(X)
    plot_correlation_heatmap(corr_matrix)
    
    # Calculate and print VIF
    vif_results = calculate_vif(X)
    print("Variance Inflation Factors:")
    print(vif_results)
    
    # Identify highly correlated features
    high_correlation = np.abs(corr_matrix) > 0.8
    high_corr_pairs = [(corr_matrix.index[i], corr_matrix.columns[j]) 
                       for i in range(len(corr_matrix.index)) 
                       for j in range(i+1, len(corr_matrix.columns)) 
                       if high_correlation.iloc[i,j] and i != j]
    
    print("\nHighly correlated feature pairs:")
    for pair in high_corr_pairs:
        print(f"{pair[0]} and {pair[1]}: {corr_matrix.loc[pair[0], pair[1]]:.2f}")

if __name__ == "__main__":
    main()

