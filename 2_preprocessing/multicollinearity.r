# Load necessary libraries
library(ggplot2)
library(ggcorrplot)
library(car) # For VIF
library(dplyr)
library(readr) # For reading CSV files

source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")


# Function to calculate correlation matrix
calculate_correlation_matrix <- function(df) {
    return(cor(df, use = "pairwise.complete.obs"))
}

# Function to plot correlation heatmap
plot_correlation_heatmap <- function(corr_matrix) {
    p <- ggcorrplot(corr_matrix,
        hc.order = TRUE,
        type = "lower",
        lab = TRUE,
        lab_size = 4, # Increase label size
        colors = c("blue", "white", "red"), # Adjust color scheme
        title = "Correlation Heatmap",
        ggtheme = ggplot2::theme_minimal(base_size = 14), # Increase base font size
        tl.cex = 1.2 # Increase text size for variable names
    )

    # Save the plot as a high-resolution JPEG
    jpeg("0_results/correlation_heatmap.jpg", width = 1500, height = 1000, quality = 100)
    print(p)
    dev.off()

    # Print the plot in a separate window (for RStudio or VSCode)
    print(p)
}

# Function to calculate Variance Inflation Factor (VIF)
calculate_vif <- function(df) {
    # Standardize the dataset
    df_scaled <- as.data.frame(scale(df))

    # Fit a linear model with all predictors
    vif_values <- vif(lm(1:nrow(df_scaled) ~ ., data = df_scaled))

    # Convert to data frame
    vif_df <- data.frame(Feature = names(vif_values), VIF = vif_values)
    vif_df <- vif_df[order(-vif_df$VIF), ]

    return(vif_df)
}

# Function to identify highly correlated features
find_highly_correlated <- function(corr_matrix, threshold = 0.8) {
    high_corr_pairs <- which(abs(corr_matrix) > threshold, arr.ind = TRUE)
    high_corr_pairs <- high_corr_pairs[high_corr_pairs[, 1] != high_corr_pairs[, 2], ] # Remove self-correlations

    unique_pairs <- unique(t(apply(high_corr_pairs, 1, sort))) # Ensure pairs are unique

    if (nrow(unique_pairs) == 0) {
        message("No highly correlated feature pairs found.")
    } else {
        message("Highly correlated feature pairs:")
        for (i in 1:nrow(unique_pairs)) {
            cat(
                rownames(corr_matrix)[unique_pairs[i, 1]], "and", colnames(corr_matrix)[unique_pairs[i, 2]],
                ": ", round(corr_matrix[unique_pairs[i, 1], unique_pairs[i, 2]], 2), "\n"
            )
        }
    }
}

# Main function to run the data preparation steps

# Load and preprocess the dataset (modify the path as needed)
secondary_CMIP6 <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 10000, asymptote = "quarter_biomass")

terraclim_pars <- c("mean_srad", "mean_soil", "mean_temp", "mean_vpd", "mean_aet", "mean_def", "mean_pdsi")
cmip6_pars <- c("nssh", "musc", "sdsr", "nsat")

soil <- c("nitro", "phh2o", "ocd", "cec", "sand", "clay", "soc", "cfvo")


df <- secondary_CMIP6 %>%
    select(c("age", "asymptote", "biomass", soil, terraclim_pars))
 
# Calculate correlation matrix
corr_matrix <- calculate_correlation_matrix(df)

# Plot correlation heatmap
plot_correlation_heatmap(corr_matrix)

# Compute VIF values
vif_results <- calculate_vif(df)
print("Variance Inflation Factors:")
print(vif_results)

# Identify highly correlated features
find_highly_correlated(corr_matrix, threshold = 0.3)
colnames(secondary_CMIP6)
