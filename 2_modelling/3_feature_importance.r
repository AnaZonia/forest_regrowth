# Load required libraries
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

# Source other scripts
source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_feature_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------ Calculate Permutation Importance ---------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

calculate_permutation_importance <- function(model, data, data_pars) {
    data_pars <- c(unlist(data_pars), "age")
    # Get baseline performance
    baseline_pred <- growth_curve(model$par, data, model$par["lag"])
    baseline_r2 <- calc_r2(data, baseline_pred)

    # Identify categorical variables
    categorical_vars <- list(
        ecoregion = grep("^ecoreg_", data_pars, value = TRUE),
        topography = grep("^topography_", data_pars, value = TRUE)
    )

    # Calculate importance for each variable
    importance <- numeric(length(data_pars))
    names(importance) <- data_pars

    for (i in seq_along(data_pars)) {
        var_name <- data_pars[i]

        # Skip if already permuted as part of a categorical variable
        if (var_name %in% unlist(categorical_vars)) next

        # Create permuted dataset
        data_permuted <- data

        if (var_name %in% unlist(categorical_vars$ecoregion)) {
            # Permute all ecoregion dummy variables together
            ecoregion_vars <- categorical_vars$ecoregion
            permuted_indices <- sample(nrow(data))
            data_permuted[ecoregion_vars] <- data[ecoregion_vars][permuted_indices, ]
        } else if (var_name %in% unlist(categorical_vars$topography)) {
            # Permute all topography dummy variables together
            topography_vars <- categorical_vars$topography
            permuted_indices <- sample(nrow(data))
            data_permuted[topography_vars] <- data[topography_vars][permuted_indices, ]
        } else {
            # Permute single variable
            data_permuted[[var_name]] <- sample(data[[var_name]])
        }

        # Get performance with permuted variable
        permuted_pred <- growth_curve(model$par, data_permuted, model$par["lag"])
        permuted_r2 <- calc_r2(data_permuted, permuted_pred)

        # Importance = decrease in performance
        importance[i] <- baseline_r2 - permuted_r2
    }

    # Create dataframe for plotting
    importance_df <- data.frame(
        variable = names(importance),
        importance = importance,
        importance_pct = 100 * importance / sum(importance)
    )

    # Sort by importance
    importance_df <- importance_df[order(-importance_df$importance), ]

    return(importance_df)
}

# ----------- Main Analysis ---------------------

all_results <- list()
groups <- c("full_amazon", "nearest_mature")
group_names <- c("Amazon-wide Average", "Nearest Neighbor")

for (i in seq_along(groups)) {
    asymptote <- groups[i]
    data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000, asymptote = asymptote)
    norm_data <- normalize_independently(data)$train_data
    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
    r2 <- calc_r2(norm_data, pred)
    importance_results <- calculate_permutation_importance(model, norm_data, data_pars)
    importance_results$importance_scaled <- importance_results$importance_pct * r2 / 100
    importance_results$group <- group_names[i]
    all_results[[i]] <- importance_results
}

# Combine importance results
all_data <- bind_rows(all_results) %>%
    filter(importance_pct > 1)

all_data$group <- factor(
    all_data$group,
    levels = group_names
)

# Map variable short names to full names
variable_names <- c(
    age = "Age",
    sur_cover = "Surrounding Mature Forest Cover",
    num_fires = "Number of Fires",
    dist = "Distance to Nearest Mature Forest",
    indig = "Indigenous Area",
    protec = "Protected Area",
    floodable_forests = "Floodable Forests",
    phh2o = "Soil pH",
    sand = "Sand Content",
    ocs = "Organic Carbon Stock",
    ocd = "Organic Carbon Density",
    cfvo = "Coarse Fragments Volume",
    nitro = "Soil Nitrogen",
    cec = "Cation Exchange Capacity",
    clay = "Clay Content",
    soc = "Soil Organic Carbon",
    mean_pr = "Mean Precipitation",
    mean_temp = "Mean Temperature",
    mean_srad = "Mean Solar Radiation",
    mean_def = "Mean Climatic Water Deficit",
    mean_vpd = "Mean Vapor Pressure Deficit",
    mean_aet = "Mean Actual Evapotranspiration",
    mean_soil = "Mean Soil Moisture",
    mean_pdsi = "Mean Palmer Drought Severity Index"
)
all_data$variable <- factor(all_data$variable, levels = names(variable_names))

# Custom colors as before
custom_colors <- c(
    "#003f5c", "#2f4b7c", "#665191", "#a05195",
    "#d45087", "#f95d6a", "#ff7c43", "#ffa600", "#ffc300", "#ffda6a"
)
n_vars <- length(unique(all_data$variable))
recycled_colors <- rep(custom_colors, length.out = n_vars)

# Plot (final in-memory results)
p <- ggplot(all_data, aes(x = group, y = importance_scaled, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(
        labels = group_names,
        name = "Asymptote"
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.05)),
        name = "R²"
    ) +
    scale_fill_manual(
        values = recycled_colors,
        labels = variable_names,
        name = "Variable"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "right",
        text = element_text(size = 12),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )

ggsave(
    "./0_results/figures/model_performance.jpg",
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
)
