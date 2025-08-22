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



theta_list <- list()
lag_list <- list()

for (i in range(1:5)) {
        
    set.seed(i)

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

        
    # Compute age frequencies
    age_freq <- field_data %>%
        count(age, name = "freq")

    # Add frequency to df
    field_data <- field_data %>%
        left_join(age_freq, by = "age") %>%
        mutate(weight = 1 / freq)

    # For each site, sample one measurement, with higher chances for rarer ages
    field_data <- field_data %>%
        group_by(site_id) %>%
        slice_sample(n = 1, weight_by = weight) %>%
        ungroup() %>%
        select(-biome, -lat, -lon, -site_id, -plot_id) %>%
            mutate(across(any_of(categorical), as.factor))

    field_data <- dummy_cols(field_data,
        select_columns = categorical,
        remove_first_dummy = TRUE,
        remove_selected_columns = TRUE
    )
    # remove points with biomass > 400
    field_data <- field_data %>%
        filter(biomass <= 200)


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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # --------------- Get theta from field data --------------- #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # - We use only age to find the shape parameter
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


    field_data_scaled <- apply_min_max_scaling(field_data, train_stats)
    common_cols <- intersect(colnames(field_data_scaled), colnames(sat_data_scaled))
    field_data_scaled <- field_data_scaled[, common_cols]
    sat_data_scaled <- sat_data_scaled[, common_cols]
    sat_data_scaled$satellite <- 1
    field_data_scaled$satellite <- 0
    unite_field_satellite <- rbind(field_data_scaled, sat_data_scaled)

    init_params <- find_combination_pars(
        basic_pars = c("theta", "lag", "k0"),
        data_pars = data_pars_options(colnames(unite_field_satellite))[["all_mean_climate"]],
        unite_field_satellite
    )

    model <- run_optim(unite_field_satellite, init_params, conditions)

    theta_list[[i]] <- model[["par"]]["theta"]
    lag_list[[i]] <- model[["par"]]["lag"]

}

# run this 5 times and get mean and sd of theta and lag



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Get R2 of sat_model when predicting field data ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



init_params <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(satellite_data))[["all_mean_climate"]],
    sat_data_scaled
)

sat_fit_result <- run_optim(sat_data_scaled, init_params, conditions)
sat_pred_biomass <- growth_curve(sat_fit_result$par, data = sat_data_scaled, lag = sat_fit_result$par["lag"])
r2_satellite <- calc_r2(sat_data_scaled, sat_pred_biomass)
r2_satellite


field_pred_biomass <- growth_curve(field_fit_result$par, data = field_data_scaled)

r2_field <- calc_r2(field_data_scaled, field_pred_biomass)
r2_field


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Exporting results ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# CSV with the R2 and the theta value
write.csv(data.frame(
    r2_field = r2_field,
    r2_satellite = r2_satellite,
    theta = theta
), file = "./0_results/0_field_results.csv", row.names = FALSE)





# Assuming pred = predicted AGB, obs = norm_data$biomass
df <- data.frame(
    Predicted = field_pred_biomass,
    Observed = field_data_scaled$biomass
)

ext <- ggplot(df, aes(x = Predicted, y = Observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 2) +
    labs(
        x = "Predicted Biomass (Mg/ha)",
        y = "Observed Biomass (Mg/ha)"
    ) +
    coord_cartesian(expand = FALSE) +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black", size = 28, family = "Helvetica"),
        axis.text = element_text(color = "black", size = 18, family = "Helvetica"),
        legend.position = "none"
    )

# Save to file
ggsave("./0_results/figures/extended/predicted_vs_observed_field.png",
    plot = ext
)






# Histogram of field ages

ext <- ggplot(field_data, aes(x = age)) +
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

# Save to file
ggsave("./0_results/figures/extended/field_age_histogram.png",
    plot = ext, width = 1800, height = 1400, res = 300
)
