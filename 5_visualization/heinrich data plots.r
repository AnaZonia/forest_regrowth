library(tidyverse)
library(ggplot2)

# Load Data
df <- read.csv("~/Documents/data/mapbiomas_heinrich_lulc_regions.csv")
df <- subset(df, amazon_quarter_regions == 1)

# Compute mean biomass per age group
mean_biomass <- df %>%
    group_by(mapbiomas_age) %>%
    summarise(
        biomass = mean(biomass, na.rm = TRUE),
        sd_biomass = sd(biomass, na.rm = TRUE)
    )

# Normalize biomass (optional)
# mean_biomass$biomass <- mean_biomass$biomass - min(mean_biomass$biomass)

# Define Chapman-Richards function
chapman_richards <- function(age, theta2, theta3) {
    271 * (1 - exp(-theta2 * age))^theta3 + min(mean_biomass$biomass)
}

# Fit the model using Nonlinear Least Squares
fit <- nls(
    biomass ~ chapman_richards(mapbiomas_age, theta2, theta3),
    data = mean_biomass,
    start = list(theta2 = 0.001, theta3 = 0.5) # Initial guesses
)

# Extract fitted parameters
theta2_fit <- coef(fit)["theta2"]
theta3_fit <- coef(fit)["theta3"]

# Generate predictions using fitted parameters
mean_biomass$predicted_biomass <- chapman_richards(mean_biomass$mapbiomas_age, theta2_fit, theta3_fit)

# Sample 15 random biomass values per age
sampled_biomass <- df %>%
    group_by(mapbiomas_age) %>%
    sample_n(size = min(15, n()), replace = FALSE) %>% # Ensure we don't sample more than available
    ungroup()

# Plot
ggplot() +
    # Scatter plot of sampled biomass values (light blue)
    geom_jitter(
        data = sampled_biomass, aes(x = mapbiomas_age, y = biomass),
        color = "#65bedb", alpha = 0.5, size = 2, width = 0.3
    ) +
    # Scatter plot of means (dark blue)
    geom_point(
        data = mean_biomass, aes(x = mapbiomas_age, y = biomass),
        size = 3, color = "blue", alpha = 0.8
    ) +
    # SD error bars
    geom_errorbar(
        data = mean_biomass, aes(x = mapbiomas_age, ymin = biomass - sd_biomass, ymax = biomass + sd_biomass),
        width = 1, color = "blue", alpha = 0.5
    ) +
    # Fitted curve (red)
    geom_line(
        data = mean_biomass, aes(x = mapbiomas_age, y = predicted_biomass),
        color = "red", size = 1.5
    ) +
    # Labels
    labs(
        title = "Chapman-Richards Growth Model Fit",
        x = "Forest Age (years)",
        y = "Biomass (Mg/ha)"
    ) +
    theme_minimal(base_size = 20) + # Larger font for readability
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"), # Centered title
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.text = element_text(size = 16)
    )
