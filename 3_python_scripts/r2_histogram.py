import pandas as pd


df = pd.read_csv("~/Documents/forest_regrowth/0_results/all_iters.csv")
# df1 = pd.read_csv("~/Documents/forest_regrowth/0_results/non_all_export.csv")
# df3 = pd.read_csv("~/Documents/forest_regrowth/0_results/aggregated_all_1985_forest.csv")
# df4 = pd.read_csv("~/Documents/forest_regrowth/0_results/keep_only_land_use_and_landscape.csv")
# df5 = pd.read_csv("~/Documents/forest_regrowth/0_results/biome_interval_aggregated.csv")
# df = pd.concat([df, df1, df3, df4, df5])
# Calculate mean r2_final per each value of datasource
mean_r2_datasource = df.groupby('datasource')['r2_final'].mean().reset_index()

# Calculate mean r2_final per each value of biome
mean_r2_biome = df.groupby('biome')['r2_final'].mean().reset_index()

# Calculate mean r2_final per each value of remove_land_use
mean_r2_remove_land_use = df.groupby('land_use')['r2_final'].mean().reset_index()

# Calculate mean r2_final per each value of remove_landscape
mean_r2_remove_landscape = df.groupby('landscape')['r2_final'].mean().reset_index()

# Display results
print(mean_r2_datasource)
print(mean_r2_biome)
print(mean_r2_remove_land_use)
print(mean_r2_remove_landscape)
