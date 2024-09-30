# main.py
from run import load_and_preprocess_data
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

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create a color palette
    palette = sns.color_palette("viridis", n_colors=len(mean_values_by_biome))
    color_dict = dict(zip(mean_values_by_biome.index, palette))

    sns.histplot(data=plot_data, x='value', hue='biome', multiple='stack', kde=True, palette=color_dict, ax=ax)

    ax.set_title(title)
    ax.set_xlabel('Biomass (tons per hectare)')
    ax.set_ylabel('Frequency')
    
    # Create custom legend
    legend_elements = [Line2D([0], [0], color=color_dict[biome], lw=4, 
                              label=f'{biome} (Mean: {mean:.2f})') 
                       for biome, mean in mean_values_by_biome.items()]
    
    ax.legend(handles=legend_elements, title='Biome', title_fontsize='13', fontsize='10', 
              loc='center left', bbox_to_anchor=(1, 0.5))

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

    fig, ax = plt.subplots(figsize=(12, 8))
    
    palette = sns.color_palette("viridis", n_colors=len(df['biome'].unique()))
    
    for i, biome in enumerate(df['biome'].unique()):
        biome_data = biomass_stats[biomass_stats['biome'] == biome]
        
        # Plot mean line
        ax.plot(biome_data['age'], biome_data['mean'], 
                label=biome, color=palette[i], linewidth=2)
        
        # Add semi-transparent buffer area
        ax.fill_between(biome_data['age'], 
                        biome_data['mean'] - biome_data['std'],
                        biome_data['mean'] + biome_data['std'],
                        color=palette[i], alpha=0.3)

    ax.set_title('Mean Biomass by Age for Each Biome (with Standard Deviation)', fontsize=16)
    ax.set_xlabel('Age (years)', fontsize=12)
    ax.set_ylabel('Mean Biomass (tons per hectare)', fontsize=12)
    
    ax.legend(title='Biome', title_fontsize='13', fontsize='10', loc='center left', bbox_to_anchor=(1, 0.5))
    ax.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    return fig


def plot_amazon_quartile_comparison(data, variable='nearest_mature'):
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
    fig, ax = plt.subplots(figsize=(12, 8))

    # Colors for each quartile
    colors = ['blue', 'green', 'orange', 'red']

    # Plot for each quartile
    for i, (stats, label) in enumerate(zip([q1_stats, q4_stats], 
                                           ['1st Quartile', '4th Quartile'])):
        ax.plot(stats['age'], stats['mean'], label=label, color=colors[i])
        ax.fill_between(stats['age'], 
                        stats['mean'] - stats['std'],
                        stats['mean'] + stats['std'],
                        color=colors[i], alpha=0.3)

    ax.set_title(f'Mean Biomass by Age for Amazon: {variable.capitalize()} Quartile Comparison', fontsize=16)
    ax.set_xlabel('Age (years)', fontsize=12)
    ax.set_ylabel('Mean Biomass (tons per hectare)', fontsize=12)
    
    ax.legend(title=f'{variable.capitalize()} Quartile', title_fontsize='13', fontsize='10', loc='upper left')
    ax.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    return fig


def main():
    pars = [
        "nearest_mature", "age", "lulc_sum_21", "lulc_sum_15", "lulc_sum_39",
        "lulc_sum_41", "num_fires_before_regrowth", "sur_cover",
        "cwd"
    ]

    X, y, _, _, _ = load_and_preprocess_data("./0_data/non_aggregated_100k_all.csv", pars, keep_all_data=True)
    
    data = {'X': X, 'y': y}

    # biome_names = {
    #     1: 'Amazon',
    #     4: 'Atlantic'
    # }

    # # Create histograms
    # hist_1 = create_histogram(data, age=1, biome_names=biome_names, 
    #                           title='Distribution of biomass by Biome for 1 year old forests')
    # hist_30 = create_histogram(data, age=30, biome_names=biome_names, 
    #                            title='Distribution of biomass by Biome for 30 year old forests')

    # # Create mean plot with standard deviation
    # mean_std_plot = create_mean_std_plot(data, biome_names)

    # # Create Amazon nearest_mature comparison plot
    # mature = plot_amazon_quartile_comparison(data, biome_names, variable='nearest_mature')
    cwd = plot_amazon_quartile_comparison(data, variable='lulc_sum_39')

    plt.show()

if __name__ == "__main__":
    main()
