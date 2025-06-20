

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(cowplot) # For legend extraction
library(ggpubr) # For legend extraction

# Set global options
options(stringsAsFactors = FALSE)
theme_set(theme_minimal(base_size = 20))

# ---------------------------- Data Loading ----------------------------
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

lag_data <- read.csv("0_results/pred_vs_obs_amazon_lag.csv")
intercept_data <- read.csv("0_results/pred_vs_obs_amazon_intercept.csv")

# Calculate the age lag
model <- readRDS("0_results/amazon_model_lag.rds")
lag <- round(model$par["lag"])

# ---------------------------- Data Preparation ----------------------------

intercept_summary <- intercept_data %>%
    group_by(age) %>%
    summarise(
        mean_pred_intercept = median(pred, na.rm = TRUE),
        sd_pred_intercept = sd(pred, na.rm = TRUE)
    )

lag_summary <- lag_data %>%
    group_by(age) %>%
    summarise(
        mean_pred_lag = median(pred, na.rm = TRUE),
        sd_pred_lag = sd(pred, na.rm = TRUE)
    )

# in both, keep only sd_pred_lag for ages greater than lag
intercept_summary <- intercept_summary %>%
    mutate(sd_pred_intercept = if_else(age < lag+1, 0, sd_pred_intercept))

lag_summary <- lag_summary %>%
    mutate(sd_pred_lag = if_else(age < lag+1, 0, sd_pred_lag))

satellite_summary <- lag_data %>%
    group_by(obs_age) %>%
    summarise(
        mean_obs = median(obs, na.rm = TRUE),
        sd_obs = sd(obs, na.rm = TRUE)
    ) %>%
    rename(age = obs_age)

# First, combine all data including future predictions
all_pred_data <- full_join(aggregated_field, lag_summary, by = "age") %>%
    full_join(intercept_summary, by = "age") %>%
    full_join(satellite_summary, by = "age")

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

    # Field data points
    geom_point(
        data = aggregated_field,
        aes(x = age, y = field_biomass, color = "Field Measurements"),
        size = 3, alpha = 0.7
    ) +

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
        "text", x = (lag + 1) + 2, y = 320,
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

