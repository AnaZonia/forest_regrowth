# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                 Tests for multicollinearity
#
#                  Ana Avila - August 2025
#
#  Tests for multicollinearity with VIF.
#  Multicollinear variables are subsequently removed in 1_parameters.r
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(ggcorrplot)
library(car) # For VIF
library(dplyr)
library(readr)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------- Functions -------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

calculate_correlation_matrix <- function(df) {
    return(cor(df, use = "pairwise.complete.obs"))
}


plot_correlation_heatmap <- function(corr_matrix) {
    p <- ggcorrplot(corr_matrix,
        hc.order = TRUE,
        type = "lower",
        lab = TRUE,
        lab_size = 4,
        colors = c("blue", "white", "red"),
        title = "Correlation Heatmap",
        ggtheme = ggplot2::theme_minimal(base_size = 14),
        tl.cex = 1.2
    )

    print(p)
}

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------- Verify multicollinearity ------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

csv_files <- list.files("./0_data/grid_10k_amazon_secondary", pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(read_csv) %>%
    bind_rows() %>%
    filter(biome == 1) %>%
    select(-c("lon", "lat", "quarter_biomass", "ecoreg_biomass", "quarter", "secondary_area", "biome")) %>%
    select(-c(
        # "mean_def"
    ))
    # List to update variables found to be multicollinear

predictors <- setdiff(colnames(df), "biomass")

mod <- lm(as.formula(paste("biomass ~", paste(predictors, collapse = " + "))), data = df)

vif_df <- data.frame(Feature = names(vif(mod)), VIF = vif(mod))
vif_df <- vif_df[order(-vif_df$VIF), ]
print(vif_df)

corr_matrix <- calculate_correlation_matrix(df)

# Identify highly correlated features
find_highly_correlated(corr_matrix, threshold = 0.4)
