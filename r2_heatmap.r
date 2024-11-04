library(tidyverse)
library(ggplot2)

# Load and preprocess data
# data <- read.csv("./new_data_yearly/non_aggregated_final_results.csv")
data <- read.csv("./new_data_yearly/aggregated_direct_results.csv")

# Convert interval and biome to factors with specified levels for custom ordering
data$interval <- factor(data$interval, levels = c("all", "5yr", "10yr", "15yr"))
data$biome <- factor(data$biome, levels = c("amaz", "atla", "both"))
data$model <- ifelse(data$model == "optim", data$basic_pars, data$model)

# Create a combined column for biome and interval
data <- data %>%
    mutate(biome_interval = paste(biome, interval, sep = " - "))

# Define a colorblind-friendly palette for heatmap colors
heat_colors <- hcl.colors(12, palette = "Heat")

# Plot the heatmap with improved labels and centralized title
ggplot(data, aes(x = biome_interval, y = model, fill = r2_final)) +
    geom_tile(color = "white") + # White borders for better tile separation
    geom_text(aes(label = round(r2_final, 2)), color = "black", size = 3) + # Display rounded R^2 values
    scale_fill_gradientn(colors = heat_colors, name = expression(R^2)) + # Color scale with R^2 symbol

    # Updated labels and title
    labs(
        title = "Model Performance Across Biomes and Temporal Intervals",
        subtitle = "R-Squared Values (R²) of Models for Different Biomes and Aggregation Periods",
        x = "Biome and Land Use Interval",
        y = "Model Type"
    ) +

    # Aspect ratio and theme for clarity
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Adjust x-axis text for readability
        axis.text.y = element_text(size = 12), # Larger font size for y-axis
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "lines") # Adjust legend key size for readability
    )


# Define a function to calculate mean R-squared for a specified grouping variable
calculate_mean_r2 <- function(group_var) {
    data %>%
        group_by(.data[[group_var]]) %>%
        summarise(mean_r2_final = mean(r2_final, na.rm = TRUE)) %>%
        rename(Group = .data[[group_var]]) %>%
        ungroup()
}

# Apply the function to each grouping variable and print results
for (group in c("biome", "model", "data")) {
    mean_r2 <- calculate_mean_r2(group)
    cat("Mean R-squared per", group, ":\n")
    print(mean_r2)
    cat("\n")
}
