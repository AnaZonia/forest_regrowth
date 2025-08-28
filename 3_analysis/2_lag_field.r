# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#       Plot the lag-corrected and uncorrected models
#
#                 Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(foreach)
library(doParallel)
library(tidyverse)
library(ggplot2)
library(cowplot) # For legend extraction
library(ggpubr) # For legend extraction

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)
# For plotting
options(stringsAsFactors = FALSE)
theme_set(theme_minimal(base_size = 20))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Model fitting and prediction
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)
norm_data <- normalize_independently(data)
norm_data <- norm_data$train_data

predictions <- data.frame(age = 1:200)

for (basic_pars_name in names(basic_pars_options)) {
    basic_pars <- basic_pars_options[[basic_pars_name]]
    norm_data_iter <- norm_data

    if (basic_pars_name == "intercept") {
        # Force the data to intercept through zero
        mean_biomass_at_zero_age <- median(norm_data_iter$biomass[norm_data_iter$age == 1], na.rm = TRUE)
        norm_data_iter$biomass <- norm_data_iter$biomass - mean_biomass_at_zero_age
    }

    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data_iter)
    model <- run_optim(norm_data_iter, init_pars, conditions)

    if (basic_pars_name == "lag") {
        lag <- round(model$par["lag"])
    }

    biomass_df <- data.frame(matrix(nrow = nrow(norm_data_iter), ncol = 0))
    for (age in 1:nrow(predictions)) {
        norm_data_iter$age <- age
        pred <- growth_curve(model$par, norm_data_iter)
        biomass_df[[as.character(age)]] <- pred
    }

    mean_biomass <- apply(biomass_df, 2, mean, na.rm = TRUE)
    sd_biomass <- apply(biomass_df, 2, sd, na.rm = TRUE)
    
    predictions[[paste0("mean_", basic_pars_name)]] <- mean_biomass
    predictions[[paste0("sd_", basic_pars_name)]] <- sd_biomass
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load and summarize field and satellite data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
field_data <- read.csv("0_data/groa_field/field_predictors.csv")
field_data <- subset(field_data, biome == 1) # Filter for Amazon biome
field_data <- field_data %>%
    rename(biomass = field_biom) %>%
    mutate(age = floor(age + 0.5))

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
    ungroup()

# # get average satellite biomass per age
# field_data <- field_data %>%
#     select(age, biomass) %>%
#     group_by(age) %>%
#     summarise(
#         mean_obs = mean(biomass, na.rm = TRUE),
#         sd_obs = sd(biomass, na.rm = TRUE)
#     ) %>%
#     mutate(age = age + lag)



# get average satellite biomass per age
aggregated_satellite <- norm_data %>%
    select(age, biomass) %>%
    group_by(age) %>%
    summarise(
        mean_obs = mean(biomass, na.rm = TRUE),
        sd_obs = sd(biomass, na.rm = TRUE)
    ) %>%
    mutate(age = age + lag)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


pred_plot <- subset(predictions, age >= 0 & age <= 75) #max(aggregated_satellite$age)

# For pred_plot, clamp the ribbon lower bound at 0
pred_plot <- pred_plot %>%
    mutate(
        ymin_lag = pmax(mean_lag - sd_lag, 0),
        ymax_lag = mean_lag + sd_lag,
        ymin_intercept = pmax(mean_intercept - sd_intercept, 0),
        ymax_intercept = mean_intercept + sd_intercept
    )

# Update plot colors and linetypes to include future predictions
plot_colors <- c(
    "Lag-corrected" = "#2f4b7c",
    "Uncorrected" = "#ffa600",
    "Observed (Remote Sensing)" = "#b91523",
    "Field Measurements" = "black"
)

linetypes <- c(
    "Lag-corrected" = "solid",
    "Uncorrected" = "solid",
    "Observed (Remote Sensing)" = "21"
)


p <- ggplot() +

    # Predicted data - current
    geom_line(
        data = pred_plot,
        aes(
            x = age, y = mean_lag,
            color = "Lag-corrected",
            linetype = "Lag-corrected"
        ),
        linewidth = 1.5
    ) +

    geom_ribbon(
        data = pred_plot,
        aes(
            x = age,
            ymin = ymin_lag,
            ymax = ymax_lag,
            fill = "Lag-corrected"
        ),
        alpha = 0.2, color = NA
    ) +
    
    geom_line(
        data = pred_plot,
        aes(
            x = age,
            y = mean_intercept,
            color = "Uncorrected",
            linetype = "Uncorrected"
        ),
        linewidth = 1.5) +

    geom_ribbon(
        data = pred_plot,
        aes(
            x = age,
            ymin = ymin_intercept,
            ymax = ymax_intercept,
            fill = "Uncorrected"
        ),
        alpha = 0.2, color = NA
    ) +


    # Field data points
    geom_point(
        data = field_data,
        aes(
            x = age, y = biomass,
            color = "Field Measurements"
        ),
        size = 3, alpha = 0.7
    ) +

    # Remote sensing data
    geom_line(
        data = aggregated_satellite,
        aes(
            x = age,
            y = mean_obs,
            color = "Observed (Remote Sensing)",
            linetype = "Observed (Remote Sensing)"
        ),
        linewidth = 1.5
    )  +

    geom_ribbon(
        data = aggregated_satellite,
        aes(
            x = age,
            ymin = mean_obs - sd_obs,
            ymax = mean_obs + sd_obs,
            fill = "Observed (Remote Sensing)"
        ),
        alpha = 0.2, color = NA
    ) +



    # Vertical line for age lag
    geom_vline(
        xintercept = (lag + 1),
        linetype = "dotted", color = "black", linewidth = 1
    ) +
    annotate(
        "text",
        x = (lag + 1) + 2, y = 320,
        label = paste(lag, "year lag"),
        color = "black", size = 7, hjust = 0
    ) +

    # Scale definitions
    scale_color_manual(values = plot_colors, name = NULL) +
    scale_fill_manual(values = plot_colors, name = NULL) +
    scale_linetype_manual(values = linetypes, name = NULL) +
    # scale_y_continuous(limits = c(0, 320), expand = expansion(0,0)) +
    coord_cartesian(ylim = c(0, 330), expand = FALSE) +

    # Labels and theme
    labs(
        x = "Forest Age (years)",
        y = "Biomass (Mg/ha)"
    ) +
    theme(
        aspect.ratio = 0.5,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black", family = "Helvetica"),
        axis.text = element_text(color = "black", size = 18, family = "Helvetica"),
        legend.position = "none"
    )

p



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Save outputs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggsave(
    filename = "0_results/figures/lag_field_biomass.jpeg",
    plot = p,
    width = 15,
    height = 8,
    units = "in",
    dpi = 300
)

p <- p + theme(legend.position = "right")

# Extract legend
legend <- cowplot::get_legend(p)

# Create a blank plot with just the legend
legend_plot <- cowplot::ggdraw(legend)

# Save the legend
ggsave(
    filename = "0_results/figures/lag_field_biomass_legend.jpeg",
    plot = legend_plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 300
)




