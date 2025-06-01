# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(cowplot) # For legend extraction
library(ggpubr) # For legend extraction

install.packages("pak")
pak::pak("ggpubr")

# Set global options
options(stringsAsFactors = FALSE)
theme_set(theme_minimal(base_size = 20))

# ---------------------------- Data Loading ----------------------------
field_data <- read.csv("0_data/groa_field/aggregated_field_biomass.csv")
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

satellite_summary <- lag_data %>%
    group_by(obs_age) %>%
    summarise(
        mean_obs = median(obs, na.rm = TRUE),
        sd_obs = sd(obs, na.rm = TRUE)
    ) %>%
    rename(age = obs_age)

# First, combine all data including future predictions
all_pred_data <- full_join(field_data, lag_summary, by = "age") %>%
    full_join(intercept_summary, by = "age") %>%
    full_join(satellite_summary, by = "age")

# Update plot colors and linetypes to include future predictions
plot_colors <- c(
    "Predicted (lag-corrected)" = "blue",
    "Predicted (intercept)" = "green4",
    "Remote Sensing" = "red",
    "Field Measurements" = "black"
)

linetypes <- c(
    "Predicted (lag-corrected)" = "solid",
    "Predicted (intercept)" = "solid",
    "Remote Sensing" = "solid",
    "Field Measurements" = "solid"
)


# Modify the plotting function
p <- ggplot(all_pred_data, aes(x = age)) +

    # Remote sensing data
    geom_line(aes(y = mean_obs, color = "Remote Sensing"), size = 1) +
    geom_ribbon(
        aes(ymin = mean_obs - sd_obs, ymax = mean_obs + sd_obs, fill = "Remote Sensing"),
        alpha = 0.2, color = NA
    ) +

    # Predicted data - current
    geom_line(aes(
        y = mean_pred_lag, color = "Predicted (lag-corrected)",
        linetype = "Predicted (lag-corrected)"
    ), size = 1) +

    geom_ribbon(
        aes(
            ymin = mean_pred_lag - sd_pred_lag,
            ymax = mean_pred_lag + sd_pred_lag,
            fill = "Predicted (lag-corrected)"
        ),
        alpha = 0.2, color = NA
    ) +
    
    geom_line(aes(
        y = mean_pred_intercept, color = "Predicted (intercept)",
        linetype = "Predicted (intercept)"
    ), size = 1) +

    # geom_ribbon(
    #     aes(
    #         ymin = mean_pred_intercept - sd_pred_intercept,
    #         ymax = mean_pred_intercept + sd_pred_intercept,
    #         fill = "Predicted (intercept)"
    #     ),
    #     alpha = 0.2, color = NA
    # ) +

    # Field data points
    geom_point(
        data = field_data,
        aes(x = age, y = mean_biomass, color = "Field Measurements"),
        size = 3, alpha = 0.7
    ) +

    # Vertical line for age lag
    geom_vline(
        xintercept = (lag + 1),
        linetype = "dotted", color = "black", size = 1
    ) +
    annotate(
        "text", x = (lag + 1) + 2, y = 300,
        label = paste(lag, "year lag"),
        color = "black", size = 5, hjust = 0
    ) +
    
    # Scale definitions
    scale_color_manual(values = plot_colors, name = NULL) +
    scale_fill_manual(values = plot_colors, name = NULL) +
    scale_linetype_manual(values = linetypes, name = NULL) +
    scale_y_continuous(limits = c(0, 310), expand = expansion(mult = c(0, 0.05))) +

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
        axis.text = element_text(color = "black", size = 14, family = "Helvetica"),
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


# Extract legend
legend <- cowplot::get_legend(p)

# Create a blank plot with just the legend
legend_plot <- ggpubr::as_ggplot(legend)

# Save the legend
ggsave(
filename = "0_results/figures/lag_field_biomass_legend.jpeg",
plot = legend_plot,
width = 15,
height = 10,
units = "in",
dpi = 300
)

