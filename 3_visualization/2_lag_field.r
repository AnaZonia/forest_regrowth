



# Load necessary libraries
library(tidyverse)
library(ggplot2)



# --------------------------------- Plotting ---------------------------------#

aggregated_field <- read.csv("0_data/groa_field/aggregated_field_biomass.csv")

pred_vs_observed <- read.csv("0_results/pred_vs_obs_amazon_lag.csv")


# ---------------------------- Plotting ----------------------------


lag = pred_vs_observed$corrected_age[1] - pred_vs_observed$uncorrected_age[1]

# Compute mean values and standard deviations per age
# Restructure the data to have separate columns for each biomass type
uncorrected_data <- pred_vs_observed %>%
    group_by(corrected_age) %>%
    summarise(
        mean_pred_uncorrected = median(pred_uncorrected, na.rm = TRUE),
        sd_pred_uncorrected = sd(pred_uncorrected, na.rm = TRUE),
        mean_observed = median(obs, na.rm = TRUE),
        sd_observed = sd(obs, na.rm = TRUE)
    ) %>%
    rename(age = corrected_age)

lag_corrected_data <- pred_vs_observed %>%
    group_by(uncorrected_age) %>%
    summarise(
        mean_pred_corrected_age = median(pred_corrected_age, na.rm = TRUE),
        sd_pred_corrected_age = sd(pred_corrected_age, na.rm = TRUE),
    ) %>%
    rename(age = uncorrected_age) %>%
        filter(age <= (lag + 1))

# Merge the data frames based on age using full_join
pred_data <- full_join(lag_corrected_data, uncorrected_data, by = "age")

print(pred_data, n = 50)
colors <- c(
    "mean_pred_corrected_age" = "blue",
    "mean_pred_uncorrected" = "dodgerblue",
    "mean_observed" = "red",
    "scatter_points" = "black"
)

legend_labels <- c(
    "mean_pred_corrected_age" = "Predicted biomass before observations",
    "mean_pred_uncorrected" = "Predicted biomass",
    "mean_observed" = "Biomass estimated by remote sensing",
    "scatter_points" = "Biomass measured in field plots"
)

p <- ggplot(pred_data, aes(x = age)) +
    geom_line(aes(y = mean_observed, color = "mean_observed"), size = 1, na.rm = TRUE) +
    geom_line(aes(y = mean_pred_corrected_age, color = "mean_pred_corrected_age"), size = 1, linetype = "dotted", na.rm = TRUE) +
    geom_line(aes(y = mean_pred_uncorrected, color = "mean_pred_uncorrected"), size = 1, na.rm = TRUE) +
    geom_ribbon(aes(ymin = mean_observed - sd_observed, ymax = mean_observed + sd_observed, fill = "mean_observed"),
        alpha = 0.2, color = NA, na.rm = TRUE
    ) +
    geom_ribbon(aes(ymin = mean_pred_corrected_age - sd_pred_corrected_age, ymax = mean_pred_corrected_age + sd_pred_corrected_age, fill = "mean_pred_corrected_age"),
        alpha = 0.2, color = NA, na.rm = TRUE
    ) +
    geom_ribbon(aes(ymin = mean_pred_uncorrected - sd_pred_uncorrected, ymax = mean_pred_uncorrected + sd_pred_uncorrected, fill = "mean_pred_uncorrected"),
        alpha = 0.2, color = NA, na.rm = TRUE
    ) +
    geom_point(data = aggregated_field, aes(x = age, y = mean_biomass, color = "scatter_points"), size = 2, alpha = 0.7) +
    scale_color_manual(values = colors, name = "Legend", labels = legend_labels) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(
        x = "Forest Age (years)",
        y = "Biomass (Mg/ha)"
    ) +
    theme_minimal(base_size = 20) +
    theme(
        aspect.ratio = 1 / 2,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.x = element_text(color = "black", family = "Helvetica"),
        axis.title.y = element_text(color = "black", family = "Helvetica"),
        axis.text.x = element_text(color = "black", size = 14, family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 14, family = "Helvetica")
    ) +
    scale_y_continuous(limits = c(0, 310)) +
    geom_vline(xintercept = (lag + 1), linetype = "dotted", color = "black", size = 1) 

p <- p +
    coord_cartesian(clip = "off") + # Allow annotation outside plot panel
    annotate(
        "text",
        x = (lag + 1),
        y = 320, # Just above y-limit 310
        label = paste(lag + 1, "years"),
        hjust = 0.5, # Center horizontally on the vertical line
        size = 6,
        color = "black",
        family = "Helvetica"
    )

p

# save the plot
ggsave(
    filename = "0_results/figures/lag_field_biomass.jpeg",
    plot = p,
    width = 15,
    height = 5,
    units = "in",
    dpi = 300
)
