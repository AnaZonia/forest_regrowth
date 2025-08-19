# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#           Field Data Analysis and Model Validation
#
#                 Ana Avila - August 2025
#
#     Fit the model to the field data
#     Find theta (shape parameter) value from the field data
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(dplyr)
library(terra)
library(tidyverse)
library(foreach)
library(doParallel)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------- Field Data Cleaning ------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

field_data_full <- read.csv("./0_data/groa_field/field_predictors.csv")
field_data_full <- subset(field_data_full, biome == 1)

field_data_full <- field_data_full %>%
    rename(
        biomass = field_biom,
        asymptote = nearest_mature
    ) %>%
    mutate(age = floor(age + 0.5))

field_data <- field_data_full %>%
    group_by(site_id) %>%
    slice_sample(n = 1) %>%
    ungroup()

field_data <- field_data %>%
    select(-biome, -lat, -lon, -site_id, -plot_id) %>%
    mutate(across(any_of(categorical), as.factor))

field_data <- dummy_cols(field_data,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Train Model on Satellite Data -------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

satellite_data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)

# remove the columns in data of the form ecoreg_xxx that have the value xxx not present in field_non_repeats$ecoreg
remove_unused_dummies <- function(data, field_non_repeats, prefix) {
    cols <- grep(paste0("^", prefix, "_"), colnames(data), value = TRUE)
    valid <- paste0(prefix, "_", unique(field_non_repeats[[prefix]]))
    cols_to_remove <- setdiff(cols, valid)
    data %>% select(-all_of(cols_to_remove))
}

satellite_data <- satellite_data %>%
    remove_unused_dummies(field_data, "ecoreg") %>%
    remove_unused_dummies(field_data, "topography")

sat_data_scaled <- normalize_independently(satellite_data)
scaling_stats <- sat_data_scaled$train_stats
sat_data_scaled <- sat_data_scaled$train_data

init_params <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(satellite_data))[["all_mean_climate"]],
    sat_data_scaled
)

sat_fit_result <- run_optim(sat_data_scaled, init_params, conditions)
sat_pred_biomass <- growth_curve(sat_fit_result$par, data = sat_data_scaled, lag = sat_fit_result$par["lag"])
r2_satellite <- calc_r2(sat_data_scaled, sat_pred_biomass)
r2_satellite





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- Get theta from field data --------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# we remove other predictors (keep only the basic)
# to avoid local minima
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# for this we actually need the data without filtering one per site (full data)

field_data_scaled <- apply_min_max_scaling(field_data_full, train_stats)

init_params <- init_params[c("k0")]
init_params$theta <- 1
init_params$B0 <- 0

field_fit_result <- run_optim(field_data_scaled, init_params, conditions)
field_fit_result[["par"]]

field_pred_biomass <- growth_curve(field_fit_result$par, data = field_data_scaled)
r2 <- calc_r2(field_data_scaled, field_pred_biomass)
r2

# this value (1.8) is then used in 2_modelling/2_modelling.r
# as the default theta value



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Get R2 of sat_model when predicting field data ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# keep only the variables that are in the model
train_stats <- train_stats %>%
    filter(variable %in% names(sat_model$par))

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}

field_data_scaled <- apply_min_max_scaling(field_data, train_stats)

field_pred_biomass <- growth_curve(sat_fit_result$par, data = field_data_scaled)

r2 <- calc_r2(field_data_scaled, field_pred_biomass)
r2



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Identify and handle plots with repeated measurements
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# take only rows with more than 5 observations for the same plot_id
field_repeats <- field_data_full %>%
    group_by(plot_id) %>%
    filter(n() > 2) %>%
    ungroup()

field_repeats_scaled <- apply_min_max_scaling(field_repeats, train_stats)

pred_repeats <- growth_curve(sat_fit_result$par, data = field_repeats_scaled)

plot(field_repeats$age, field_repeats$biomass, xlab = "Age", ylab = "Biomass", main = "Predictions for site with repeated measurements")
points(field_repeats_scaled$age, pred_repeats, col = "red", pch = 19)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Extended Data Figures ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

png("./0_results/figures/extended/field_predictions_scatterplot.png", width = 800, height = 600)
plot(field_non_repeats$age, field_non_repeats$biomass, xlab = "Field Age", ylab = "Field Biomass", main = "Age vs Biomass")
points(field_non_repeats$age, pred_field, col = "red", pch = 19)
dev.off()


png("./0_results/figures/extended/predicted_vs_observed_field.png", width = 800, height = 600)
plot(pred_field, field_non_repeats$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
# add 1-1 line to the plot
abline(0, 1, col = "red", lty = 2)
dev.off()



png("./0_results/figures/extended/field_age_histogram.png", width = 1800, height = 1400, res = 300)

ggplot(field, aes(x = age)) +
    geom_histogram(
        binwidth = 5,
        fill = "grey30",
        color = "white",
        boundary = 0
    ) +
    labs(
        x = "Forest age (years)",
        y = "Number of plots"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 20),
        axis.title = element_text(face = "bold", size = 22),
        axis.ticks = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10)
    )

dev.off()
