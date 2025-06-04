# ------------------------------------------------- #
# Figure - Model Performance Section
# Plot the importance of each category of variables
# ------------------------------------------------- #

library(viridis) # for color palettes
library(ggplot2)
library(tidyverse)


importance_full_amazon <- read.csv("./0_results/importance_full_amazon.csv")
importance_nearest_mature <- read.csv("./0_results/importance_nearest_mature.csv")

importance_full_amazon$group <- "Amazon-wide Average"
importance_nearest_mature$group <- "Nearest Neighbor"

all_data <- bind_rows(importance_full_amazon, importance_nearest_mature) %>%
    filter(importance_scaled > 0)

# ---- Plot ----

# Create a mapping of short variable names to their full names
variable_names <- c(
    age = "Age",
    sur_cover = "Surrounding Mature Forest Cover",
    num_fires = "Number of Fires",
    dist = "Distance",
    indig = "Indigenous Area",
    protec = "Protected Area",
    floodable_forests = "Floodable Forests",
    phh2o = "Soil pH",
    sand = "Sand Content",
    ocs = "Organic Carbon Stock",
    ocd = "Organic Carbon Density",
    cfvo = "Coarse Fragments Volume",
    nitro = "Soil Nitrogen",
    cec = "Cation Exchange Capacity",
    clay = "Clay Content",
    soc = "Soil Organic Carbon",
    mean_pr = "Mean Precipitation",
    mean_temp = "Mean Temperature",
    mean_srad = "Mean Solar Radiation",
    mean_def = "Mean Climatic Water Deficit",
    mean_vpd = "Mean Vapor Pressure Deficit",
    mean_aet = "Mean Actual Evapotranspiration",
    mean_soil = "Mean Soil Moisture",
    mean_pdsi = "Mean Palmer Drought Severity Index"
)


# Map variable short names to full names in the data frame
all_data$variable <- factor(all_data$variable, levels = names(variable_names))

ggplot(all_data, aes(x = group, y = importance_scaled, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(
        labels = c(Amazon = "Amazon-wide Average", NearestMature = "Nearest Neighbor"),
        name = "Asymptote"
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.05)),
        name = "Feature Importance"
    ) +
    scale_fill_manual(
        values = viridis::viridis(length(unique(all_data$variable)), option = "D"),
        labels = variable_names,
        name = "Variable"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "right",
        text = element_text(size = 12),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )

# save plot as jpeg
ggsave(
    "./0_results/figures/model_performance.jpg",
    width = 10,
    height = 6,
    dpi = 300
)
