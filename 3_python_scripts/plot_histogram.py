# main.py
from run_ML import load_and_preprocess_data
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def create_histogram(data, age, biome_names, title):
    plot_data = pd.DataFrame({
        'value': data['y'][data['X']['age'] == age],
        'biome': data['X']['biome'][data['X']['age'] == age].map(biome_names)
    })

    mean_values_by_biome = plot_data.groupby('biome')['value'].mean()

    fig, ax = plt.subplots(figsize = (10, 6))
    
    # Create a color palette
    palette = sns.color_palette("viridis", n_colors = len(mean_values_by_biome))
    color_dict = dict(zip(mean_values_by_biome.index, palette))

    sns.histplot(data = plot_data, x = 'value', hue = 'biome', multiple = 'stack', kde = True, palette = color_dict, ax = ax)

    ax.set_title(title)
    ax.set_xlabel('Biomass (tons per hectare)')
    ax.set_ylabel('Frequency')
    
    # Create custom legend
    legend_elements = [Line2D([0], [0], color = color_dict[biome], lw = 4, 
                              label = f'{biome} (Mean: {mean:.2f})') 
                       for biome, mean in mean_values_by_biome.items()]
    
    ax.legend(handles = legend_elements, title = 'Biome', title_fontsize = '13', fontsize = '10', 
              loc = 'center left', bbox_to_anchor = (1, 0.5))

    plt.tight_layout()
    return fig

def create_mean_std_plot(data, biome_names):
    df = pd.DataFrame({
        'age': data['X']['age'],
        'biome': data['X']['biome'].map(biome_names),
        'biomass': data['y']
    })

    # Calculate mean and standard deviation
    biomass_stats = df.groupby(['age', 'biome'])['biomass'].agg(['mean', 'std']).reset_index()

    fig, ax = plt.subplots(figsize = (12, 8))
    
    palette = sns.color_palette("viridis", n_colors = len(df['biome'].unique()))
    
    for i, biome in enumerate(df['biome'].unique()):
        biome_data = biomass_stats[biomass_stats['biome'] == biome]
        
        # Plot mean line
        ax.plot(biome_data['age'], biome_data['mean'], 
                label = biome, color = palette[i], linewidth = 2)
        
        # Add semi-transparent buffer area
        ax.fill_between(biome_data['age'], 
                        biome_data['mean'] - biome_data['std'],
                        biome_data['mean'] + biome_data['std'],
                        color = palette[i], alpha = 0.3)

    ax.set_title('Mean Biomass by Age for Each Biome (with Standard Deviation)', fontsize = 16)
    ax.set_xlabel('Age (years)', fontsize = 12)
    ax.set_ylabel('Mean Biomass (tons per hectare)', fontsize = 12)
    
    ax.legend(title = 'Biome', title_fontsize = '13', fontsize = '10', loc = 'center left', bbox_to_anchor = (1, 0.5))
    ax.grid(True, linestyle = '--', alpha = 0.7)

    plt.tight_layout()
    return fig


def plot_amazon_quartile_comparison(data, variable = 'nearest_mature_biomass'):
    # Filter for Amazon biome
    amazon_data = pd.DataFrame({
        'age': data['X']['age'][data['X']['biome'] == 1],
        'biomass': data['y'][data['X']['biome'] == 1],
        variable: data['X'][variable][data['X']['biome'] == 1]
    })

    # Calculate quartiles
    quartiles = amazon_data[variable].quantile([0.25, 0.5, 0.75])

    # Create masks for each quartile
    q1_mask = amazon_data[variable] <= quartiles[0.25]
    q4_mask = amazon_data[variable] > quartiles[0.25]

    # Function to calculate mean and std for a group
    def calc_mean_std(group):
        return pd.Series({
            'mean': group['biomass'].mean(),
            'std': group['biomass'].std()
        })

    # Calculate mean and std for each age group in each quartile
    q1_stats = amazon_data[q1_mask].groupby('age').apply(calc_mean_std).reset_index()
    q4_stats = amazon_data[q4_mask].groupby('age').apply(calc_mean_std).reset_index()

    # Create the plot
    fig, ax = plt.subplots(figsize = (12, 8))

    # Colors for each quartile
    colors = ['blue', 'green', 'orange', 'red']

    # Plot for each quartile
    for i, (stats, label) in enumerate(zip([q1_stats, q4_stats], 
                                           ['1st Quartile', '4th Quartile'])):
        ax.plot(stats['age'], stats['mean'], label = label, color = colors[i])
        ax.fill_between(stats['age'], 
                        stats['mean'] - stats['std'],
                        stats['mean'] + stats['std'],
                        color = colors[i], alpha = 0.3)

    ax.set_title(f'Mean Biomass by Age for Amazon: {variable.capitalize()} Quartile Comparison', fontsize = 16)
    ax.set_xlabel('Age (years)', fontsize = 12)
    ax.set_ylabel('Mean Biomass (tons per hectare)', fontsize = 12)
    
    ax.legend(title = f'{variable.capitalize()} Quartile', title_fontsize = '13', fontsize = '10', loc = 'upper left')
    ax.grid(True, linestyle = '--', alpha = 0.7)

    plt.tight_layout()
    return fig

def plot_mean_biomass_comparison(df1, df2, label1='EU', label2='Non-Aggregated'):
    """
    Plot mean biomass per age for two dataframes to compare.
    
    Args:
        df1 (pd.DataFrame): First dataset containing 'age' and 'biomass' columns.
        df2 (pd.DataFrame): Second dataset containing 'age' and 'biomass' columns.
        label1 (str): Label for the first dataset.
        label2 (str): Label for the second dataset.
    """
    # Calculate mean and standard deviation by age for each dataframe
    stats1 = df1.groupby('age')['biomass'].agg(['mean', 'std']).reset_index().rename(columns={'mean': f'mean_{label1}', 'std': f'std_{label1}'})
    stats2 = df2.groupby('age')['biomass'].agg(['mean', 'std']).reset_index().rename(columns={'mean': f'mean_{label2}', 'std': f'std_{label2}'})
    
    # Merge statistics on age
    comparison_df = pd.merge(stats1, stats2, on='age', how='inner')

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot means and standard deviations for both datasets
    ax.plot(comparison_df['age'], comparison_df[f'mean_{label1}'], label=label1, color='blue', linewidth=2)
    ax.fill_between(comparison_df['age'], 
                    comparison_df[f'mean_{label1}'] - comparison_df[f'std_{label1}'],
                    comparison_df[f'mean_{label1}'] + comparison_df[f'std_{label1}'],
                    color='blue', alpha=0.3)

    ax.plot(comparison_df['age'], comparison_df[f'mean_{label2}'], label=label2, color='green', linewidth=2)
    ax.fill_between(comparison_df['age'], 
                    comparison_df[f'mean_{label2}'] - comparison_df[f'std_{label2}'],
                    comparison_df[f'mean_{label2}'] + comparison_df[f'std_{label2}'],
                    color='green', alpha=0.3)

    # Titles and labels
    ax.set_title('Comparison of Mean Biomass by Age', fontsize=16)
    ax.set_xlabel('Age (years)', fontsize=12)
    ax.set_ylabel('Mean Biomass (tons per hectare)', fontsize=12)
    
    # Legend and grid
    ax.legend(title='Dataset', title_fontsize='13', fontsize='10', loc='upper left')
    ax.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.show()

def main():

    # df1 = pd.read_csv("./0_data/eu.csv")
    # df2 = pd.read_csv("./0_data/non_aggregated_all.csv")

    # stats1 = df1.groupby('age')['biomass'].agg(['mean', 'std']).reset_index()
    # stats2 = df2.groupby('age')['biomass'].agg(['mean', 'std']).reset_index()
    # comparison_df = stats1.merge(stats2, on='age', suffixes=('_eu', '_mapbiomas'))
    # print(comparison_df)
    # print(comparison_df['std_eu'].mean())
    # print(comparison_df['std_mapbiomas'].mean())

    # plot_mean_biomass_comparison(df1, df2, label1='EU', label2='Mapbiomas')

    # df = pd.read_csv("./0_data/eu.csv")
    df = pd.read_csv("./0_data/non_aggregated_15yr.csv")

    # Select only the specified columns
    X = df[['age', 'biome', 'num_fires_before_regrowth', 'cwd']]  # Now X contains both age and biome
    y = df['biomass']

    data = {'X': X, 'y': y}

    biome_names = {
        1: 'Amazon',
        4: 'Atlantic'
    }

    # # Create histograms
    # hist_1 = create_histogram(data, age = 1, biome_names = biome_names, 
    #                           title = 'Distribution of biomass by Biome for 1 year old forests')
    # hist_30 = create_histogram(data, age = 30, biome_names = biome_names, 
    #                            title = 'Distribution of biomass by Biome for 30 year old forests')

    # Create mean plot with standard deviation
    # mean_std_plot = create_mean_std_plot(data, biome_names)
    
    # # Create Amazon nearest_mature_biomass comparison plot
    plot = plot_amazon_quartile_comparison(data, variable = 'num_fires_before_regrowth')


    plt.show()

if __name__ == "__main__":
    main()
