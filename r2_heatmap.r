# Load necessary libraries
library(dplyr)
library(tidyr)
library(scico)
set.seed(123) # For reproducible random R^2 values

# Define values for each column
data <- read.csv("./new_data_yearly/non_aggregated_final_results.csv")

# Convert interval and biome to factors to control the layout
data$data <- factor(data$data, levels = c("all", "5yr", "10yr", "15yr"))
data$biome <- factor(data$biome_name, levels = c("amaz", "atla", "both"))
data$model <- ifelse(data$model == "optim", data$basic_pars, data$model)

# Create a combined column for biome and interval to arrange as requested
data <- data %>%
    mutate(biome_data = paste(biome, data, sep = " - "))

# Use the "Heat" palette from hcl.colors for colorblind-friendly sequential coloring
heat_colors <- hcl.colors(12, palette = "Heat") # Generate 12 colors from the Heat palette

# Plot the heatmap with separate facets for aggregation types and text labels for R2 values
ggplot(data, aes(x = biome_data, y = model, fill = r2_final)) +
    geom_tile(color = "white") + # Adds white borders between tiles for clarity
    geom_text(aes(label = round(r2_final, 2)), color = "black", size = 3) + # Add R2 values in black
    scale_fill_gradientn(colors = heat_colors, name = "R^2 Value") + # Apply the Heat palette
    labs(
        title = "Heatmap of R Squared Values by Model, Biome, and Interval",
        x = "Biome - Interval",
        y = "Model Type"
    ) +
    coord_fixed(ratio = 1) + # Set aspect ratio to ensure square tiles
    # facet_wrap(~aggregation) + # Separate plots for each aggregation type
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for readability
        axis.text.y = element_text(size = 10), # Adjust y-axis font size for models
        strip.text = element_text(size = 12, face = "bold") # Bold for facet titles
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
