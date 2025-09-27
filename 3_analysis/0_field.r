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
    mutate(age = floor(age + 0.5)) %>%
    drop_na()

# Assuming field_data already cleaned as shown earlier
field_summary <- field_data %>%
    group_by(age) %>%
    summarise(
        biomass = median(biomass, na.rm = TRUE),
        .groups = "drop"
    )


plot(field_summary$age, field_summary$biomass,
    xlab = "Age (years)",
    ylab = "Biomass (Mg/ha)",
    main = "Field Data: Biomass vs Age"
)


# site_env <- field_data %>%
#     select(-age, -biomass, -plot_id) %>%
#     distinct(site_id, .keep_all = TRUE)

# # field_summary <- field_summary %>%
#     # left_join(site_env, by = "site_id")


field_data <- dummy_cols(field_data,
    select_columns = categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

field_data <- field_data %>%
    filter(biomass < 400)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Get R2 of sat_model when predicting field data ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


init_params <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(sat_data_scaled))[["all"]],
    sat_data_scaled
)

sat_fit_result <- run_optim(sat_data_scaled, init_params, conditions)
sat_pred_biomass <- growth_curve(sat_fit_result$par, data = sat_data_scaled, lag = sat_fit_result$par["lag"])
r2_satellite <- calc_r2(sat_data_scaled, sat_pred_biomass)
r2_satellite

field_pred_biomass <- growth_curve(field_fit_result$par, data = field_data_scaled)

r2_field <- calc_r2(field_data_scaled, field_pred_biomass)
r2_field
