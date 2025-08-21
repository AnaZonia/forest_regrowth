# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#         Mature Forest Biomass by Distance to Edge
#
#                    Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyverse)
library(ggplot2)

mature_biomass_distance <- read.csv("./0_data/mature_biomass/mature_biomass_distance.csv")


mature_biomass_distance <- mature_biomass_distance[,c("distance", "mature_biomass")]

ext <- ggplot(mature_biomass_distance, aes(x = distance, y = mature_biomass)) +
    geom_point() +
    labs(
        x = "Distance to Forest Edge (m)",
        y = "Mature Forest Biomass (Mg/ha)"
    ) +
    
    # Vertical line for age lag
    geom_vline(
        xintercept = 1000,
        linetype = "22", color = "red", linewidth = 1
    ) +
    annotate(
        "text",
        x = 1000, y = 350,
        label = "1km from forest edge",
        color = "black", size = 7, hjust = 0
    ) +
    
    coord_cartesian(xlim = c(0, 10000), ylim = c(0, 400), expand = FALSE) +

    theme(
        aspect.ratio = 1.5,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black", family = "Helvetica"),
        axis.text = element_text(color = "black", size = 18, family = "Helvetica"),
        legend.position = "none"
    )

ext

# Save to file
ggsave("0_results/figures/extended/mature_forest_distance_edge.jpeg",
    plot = ext,
    width = 5, height = 10, dpi = 300
)
