



# Load necessary libraries
library(tidyverse)
library(ggplot2)



# --------------------------------- Plotting ---------------------------------#

aggregated_field <- read.csv("0_data/groa_field/aggregated_field_biomass.csv")



# ---------------------------- Plotting ----------------------------



# Compute mean values and standard deviations per age
# Restructure the data to have separate columns for each biomass type
mean_biomass_data <- data %>%
    group_by(age) %>%
    summarise(
        mean_pred = median(pred, na.rm = TRUE),
        sd_pred = sd(pred, na.rm = TRUE),
        mean_biomass = median(biomass, na.rm = TRUE),
        sd_biomass = sd(biomass, na.rm = TRUE)
    ) %>%
    mutate(age = age + lag) # Adjust age by lag

pred_data <- data %>%
    group_by(age) %>%
    summarise(
        mean_pred_lag = median(pred_lag, na.rm = TRUE),
        sd_pred_lag = sd(pred_lag, na.rm = TRUE),
    ) %>%
    filter(age <= (lag + 1))

# Merge the data frames based on age using full_join
pred_data <- full_join(pred_data, mean_biomass_data, by = "age")

print(pred_data, n = 50)

colors <- c(
    "mean_pred_lag" = "blue",
    "mean_pred" = "blue",
    "mean_biomass" = "red",
    "scatter_points" = "red" # Use black color for the scatter points in the legend
)

# Define custom legend labels
legend_labels <- c(
    "mean_pred_lag" = "Predicted biomass before observations",
    "mean_pred" = "Predicted biomass",
    "mean_biomass" = "Biomass estimated by remote sensing",
    "scatter_points" = "Biomass measured in field plots" # Custom label for scatter points
)

# Create the plot
p <- ggplot(pred_data, aes(x = age)) +
    geom_line(aes(y = mean_biomass, color = "mean_biomass"), size = 1, na.rm = TRUE) +
    geom_line(aes(y = mean_pred_lag, color = "mean_pred_lag"), size = 1, linetype = "dotted") +
    geom_line(aes(y = mean_pred, color = "mean_pred"), size = 1, na.rm = TRUE) +
    geom_ribbon(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass, fill = "mean_biomass"),
        alpha = 0.2, color = NA, na.rm = TRUE
    ) +
    geom_ribbon(aes(ymin = mean_pred_lag - sd_pred_lag, ymax = mean_pred_lag + sd_pred_lag, fill = "mean_pred_lag"),
        alpha = 0.2, color = NA
    ) +
    geom_ribbon(aes(ymin = mean_pred - sd_pred, ymax = mean_pred + sd_pred, fill = "mean_pred"),
        alpha = 0.2, color = NA, na.rm = TRUE
    ) +
    scale_color_manual(values = colors, name = "Legend", labels = legend_labels) + # Custom labels
    scale_fill_manual(values = colors, guide = "none") +
    labs(
        x = "Forest Age (years)",
        y = "Biomass (Mg/ha)"
    ) +
    geom_point(data = field_aggregated, aes(x = age, y = mean_biomass, color = "scatter_points"), size = 2, alpha = 0.7) + # Scatter points, now part of the legend
    theme_minimal(base_size = 20) +
    theme(
        legend.text = element_text(size = 16, family = "Helvetica"),
        legend.title = element_blank(), # Remove legend title
        legend.position = c(0.25, 0.75),
        legend.background = element_rect(fill = "white", color = NA), # Remove background color
        legend.key = element_blank(), # Remove the black square around legend items
        aspect.ratio = 1 / 2,
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        axis.line = element_line(color = "black"), # Keep simple x and y axis lines
        axis.title.x = element_text(color = "black", family = "Helvetica"), # Set x-axis title to black, Helvetica font
        axis.title.y = element_text(color = "black", family = "Helvetica"), # Set y-axis title to black, Helvetica font
        axis.text.x = element_text(color = "black", size = 14, family = "Helvetica"), # Set x-axis labels to black, Helvetica font
        axis.text.y = element_text(color = "black", size = 14, family = "Helvetica"), # Set y-axis labels to black, Helvetica font
    ) +
    scale_y_continuous(limits = c(0, 310)) + # Set y-axis limits
    geom_vline(xintercept = (lag + 1), linetype = "dotted", color = "black", size = 1) + # Vertical dashed line at x = 34
    annotate("text", x = (lag + 2), y = 280, label = "34 years", hjust = -0.1, size = 6, color = "black", family = "Helvetica") # Add label with annotate()

p
