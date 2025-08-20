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

all_pred_results <- list()

for (basic_pars_name in names(basic_pars_options)) {
    basic_pars <- basic_pars_options[[basic_pars_name]]

    if (basic_pars_name == "intercept") {
        # Force the data to intercept through zero
        mean_biomass_at_zero_age <- median(norm_data$biomass[norm_data$age == 1], na.rm = TRUE)
        norm_data$biomass <- norm_data$biomass - mean_biomass_at_zero_age
    }

    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)

    pred_vs_obs <- data.frame(
        age = norm_data$age,
        pred = growth_curve(model$par, norm_data)
    )

    # get the biomass predictions for the ages 1 to 105 (in 35 year steps)
    for (i in c(35, 70)) {
        norm_data_future <- norm_data
        norm_data_future$age <- norm_data_future$age + i
        pred_future <- growth_curve(model$par, norm_data_future)
        df_future <- data.frame(
            age = norm_data_future$age,
            pred = pred_future
        )
        pred_vs_obs <- rbind(pred_vs_obs, df_future)
    }

    if (basic_pars_name == "intercept") {
        pred_vs_obs$pred <- round(pred_vs_obs$pred - model$par["B0"])
    } else {
        # add column obs with the age correspondent to that in norm_data
        pred_vs_obs <- cbind(pred_vs_obs, data.frame(
            obs = norm_data$biomass,
            obs_age = round(norm_data$age + model$par["lag"])
        ))
    }

    all_pred_results[[basic_pars_name]] <- list(
        model = model,
        preds = pred_vs_obs
    )
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load and summarize field and satellite data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

field <- read.csv("0_data/groa_field/field_predictors.csv")
field <- subset(field, biome == 1) # Filter for Amazon biome
field <- field %>%
    rename(field_biomass = field_biom) %>%
    mutate(age = floor(age + 0.5))

# get the average per age
aggregated_field <- field %>%
    select(site_id, age, field_biomass) %>%
    group_by(age) %>%
    summarise(field_biomass = mean(field_biomass, na.rm = TRUE))

lag_preds <- all_pred_results[["lag"]]$preds
intercept_preds <- all_pred_results[["intercept"]]$preds
lag_model <- all_pred_results[["lag"]]$model
lag <- round(lag_model$par["lag"])


# keep only sd_pred_lag for ages greater than lag
intercept_summary <- intercept_preds %>%
    group_by(age) %>%
    summarise(
        mean_pred_intercept = median(pred, na.rm = TRUE),
        sd_pred_intercept = sd(pred, na.rm = TRUE)
    ) %>%
    mutate(sd_pred_intercept = if_else(age < lag + 1, 0, sd_pred_intercept))

lag_summary <- lag_preds %>%
    group_by(age) %>%
    summarise(
        mean_pred_lag = median(pred, na.rm = TRUE),
        sd_pred_lag = sd(pred, na.rm = TRUE)
    ) %>%
    mutate(sd_pred_lag = if_else(age < lag + 1, 0, sd_pred_lag))

satellite_summary <- intercept_preds %>%
    group_by(age) %>%
    summarise(
        mean_obs = median(pred, na.rm = TRUE),
        sd_obs = sd(pred, na.rm = TRUE)
    ) %>%
    rename(age = age)

# First, combine all data including future predictions
all_pred_data <- full_join(aggregated_field, lag_summary, by = "age") %>%
    full_join(intercept_summary, by = "age") %>%
    full_join(satellite_summary, by = "age")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Update plot colors and linetypes to include future predictions
plot_colors <- c(
    "Lag-corrected" = "#2f4b7c",
    "Uncorrected" = "#ffa600",
    "Observed (Remote Sensing)" = "#f95d6a",
    "Field Measurements" = "black"
)

linetypes <- c(
    "Lag-corrected" = "solid",
    "Uncorrected" = "solid",
    "Observed (Remote Sensing)" = "dashed"
)

# Modify the plotting function
p <- ggplot(all_pred_data, aes(x = age)) +

    # # Field data points
    # geom_point(
    #     data = aggregated_field,
    #     aes(x = age, y = field_biomass, color = "Field Measurements"),
    #     size = 3, alpha = 0.7
    # ) +

    # Remote sensing data
    geom_line(aes(y = mean_obs, color = "Observed (Remote Sensing)"), linewidth = 1.5) +
    geom_ribbon(
        aes(ymin = mean_obs - sd_obs, ymax = mean_obs + sd_obs, fill = "Observed (Remote Sensing)"),
        alpha = 0.2, color = NA
    ) +

    # Predicted data - current
    geom_line(aes(
        y = mean_pred_lag, color = "Lag-corrected",
        linetype = "Lag-corrected"
    ), linewidth = 1.5) +
    geom_ribbon(
        aes(
            ymin = mean_pred_lag - sd_pred_lag,
            ymax = mean_pred_lag + sd_pred_lag,
            fill = "Lag-corrected"
        ),
        alpha = 0.2, color = NA
    ) +
    geom_line(aes(
        y = mean_pred_intercept, color = "Uncorrected",
        linetype = "Uncorrected"
    ), linewidth = 1.5) +

    # geom_ribbon(
    #     aes(
    #         ymin = mean_pred_intercept - sd_pred_intercept,
    #         ymax = mean_pred_intercept + sd_pred_intercept,
    #         fill = "Uncorrected"
    #     ),
    #     alpha = 0.2, color = NA
    # ) +

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
    scale_y_continuous(limits = c(0, 320), expand = expansion(mult = c(0, 0.05))) +

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


satellite_summary$age <- satellite_summary$age - lag
p <- ggplot(satellite_summary, aes(x = age)) +
    # Remote sensing data
    geom_line(aes(y = mean_obs, color = "Observed (Remote Sensing)"), linewidth = 1.5) +
    geom_ribbon(
        aes(
            ymin = mean_obs - sd_obs, ymax = mean_obs + sd_obs,
            fill = "Observed (Remote Sensing)"
        ),
        alpha = 0.2, color = NA
    ) +

    # Labels and theme
    labs(
        x = "Forest Age (years)",
        y = "Mean Biomass (Mg/ha)"
    )



p


head(satellite_summary)


# ---------------------------- Save Outputs ----------------------------
# Save main plot
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
