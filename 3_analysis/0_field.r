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

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------- Field Data Cleaning ------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

field_data <- read.csv("./0_data/groa_field/field_predictors.csv")
field_data <- subset(field_data, biome == 1)

field_data <- field_data %>%
    rename(
        biomass = field_biom,
        asymptote = nearest_mature
    ) %>%
    mutate(age = floor(age + 0.5))

field_data <- field_data %>%
    group_by(site_id) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
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
train_stats <- sat_data_scaled$train_stats
sat_data_scaled <- sat_data_scaled$train_data

init_params <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(satellite_data))[["all_mean_climate"]],
    sat_data_scaled
)

sat_fit_result <- run_optim(sat_data_scaled, init_params, conditions)
sat_pred_biomass <- growth_curve(sat_fit_result$par, data = sat_data_scaled, lag = sat_fit_result$par["lag"])
r2_satellite <- calc_r2(sat_data_scaled, sat_pred_biomass)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------- Get theta from field data --------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# - We remove other predictors (keep only the basic)
# to avoid local minima
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

train_stats <- train_stats %>%
    filter(variable %in% names(sat_fit_result$par))

field_data_scaled <- apply_min_max_scaling(field_data, train_stats)

init_params <- init_params[c("k0")]
init_params$theta <- 1
init_params$B0 <- 0

field_fit_result <- run_optim(field_data_scaled, init_params, conditions)
theta <- field_fit_result[["par"]]["theta"]

# this value (1.05) is then used in 2_modelling/2_modelling.r
# as the default theta value

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Get R2 of sat_model when predicting field data ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

field_data_scaled <- apply_min_max_scaling(field_data, train_stats)

theta <- field_fit_result[["par"]]["theta"]

field_pred_biomass <- growth_curve(sat_fit_result$par, data = field_data_scaled)

r2_field <- calc_r2(field_data_scaled, field_pred_biomass)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Exporting results ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# CSV with the R2 and the theta value
write.csv(data.frame(
    r2_field = r2_field,
    r2_satellite = r2_satellite,
    theta = theta
), file = "./0_results/0_field_results.csv", row.names = FALSE)


# Predicted vs. observed scatterplot
png("./0_results/figures/extended/predicted_vs_observed_field.png", width = 800, height = 600)
plot(field_pred_biomass, field_data_scaled$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
abline(0, 1, col = "red", lty = 2)
dev.off()


# Histogram of field ages
png("./0_results/figures/extended/field_age_histogram.png", width = 1800, height = 1400, res = 300)
ggplot(field_data, aes(x = age)) +
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
