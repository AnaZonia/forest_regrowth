# ------------------------------------------------- #
# Plot the variable importance scores
# ------------------------------------------------- #
tst <- tst[tst$importance_pct > 0.5, ]
# tst <- tst[tst$variable != "age", ]
tst

# Define separate lists for each category
landscape_vars <- c("sur_cover", "dist", "age")
climate_vars <- c("mean_srad", "mean_def", "mean_vpd", "mean_aet", "mean_pr", "mean_temp")
disturbance_vars <- c("num_fires")
vegetation_vars <- c("floodable_forests")
protected_vars <- c("protec", "indig")
soil_vars <- c("phh2o", "sand", "clay", "soc", "ocs", "ocd", "cfvo", "nitro", "cec", "mean_soil")

# Create a column to store the category for each variable
tst <- tst %>%
    mutate(category = case_when(
        variable %in% landscape_vars ~ "Landscape",
        variable %in% climate_vars ~ "Climate",
        variable %in% disturbance_vars ~ "Disturbance",
        variable %in% vegetation_vars ~ "Vegetation",
        variable %in% protected_vars ~ "Protected Area",
        variable %in% soil_vars ~ "Soil",
        TRUE ~ "Other" # Fallback if variable doesn't match
    ))

# Define custom colors for each category
category_colors <- c(
    "Landscape" = "#F0E442",
    "Climate" = "#0072B2",
    "Disturbance" = "#CC79A7",
    "Vegetation" = "#009E73",
    "Protected Area" = "#E69F00",
    "Soil" = "#D55E00"
)

# Create a mapping of short variable names to their full names
variable_names <- c(
    age = "Age",
    sur_cover = "Surrounding Mature Forest Cover",
    mean_srad = "Mean Solar Radiation",
    mean_def = "Mean Climatic Water Deficit",
    num_fires = "Number of Fires",
    phh2o = "Soil pH",
    mean_vpd = "Mean Vapor Pressure Deficit",
    mean_aet = "Mean Actual Evapotranspiration",
    floodable_forests = "Floodable Forests",
    sand = "Sand Content",
    mean_soil = "Mean Soil Moisture",
    protec = "Protected Area",
    ocs = "Organic Carbon Stock",
    ocd = "Organic Carbon Density",
    cfvo = "Coarse Fragments Volume",
    nitro = "Soil Nitrogen",
    dist = "Distance",
    indig = "Indigenous Area",
    cec = "Cation Exchange Capacity",
    clay = "Clay Content",
    mean_pr = "Mean Precipitation",
    mean_temp = "Mean Temperature",
    soc = "Soil Organic Carbon"
)

# Add full names to the importance dataframe
tst$full_name <- variable_names[tst$variable]

# Create the plot with color-coded categories
importance_plot <- ggplot(tst, aes(x = reorder(full_name, importance), y = importance, fill = category)) +
    geom_col(width = 0.9) + # Reduce bar spacing
    scale_fill_manual(values = category_colors) + # Apply custom colors
    coord_flip() +
    labs(
        y = "Importance Score",
        fill = "Category"
    ) +
    theme_minimal(base_size = 16) + # Increases overall text size
    theme(
        legend.text = element_text(size = 16, family = "Helvetica"),
        legend.title = element_blank(), # Remove legend title
        legend.position = c(0.75, 0.25),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        axis.line = element_line(color = "black"), # Keep simple x and y axis lines
        axis.title.y = element_blank(), # Bigger X-axis title
        axis.text.x = element_text(color = "black", size = 14, family = "Helvetica"), # Set x-axis labels to black, Helvetica font
        axis.text.y = element_text(color = "black", size = 14, family = "Helvetica"), # Set y-axis labels to black, Helvetica font
    ) +
    scale_y_continuous(expand = c(0, 0)) # Remove space between bars and y-axis



# Print the plot
print(importance_plot)


# Optionally, save the plot to a file
ggsave("variable_importance_plot.png", plot = importance_plot, width = 10, height = 6)
