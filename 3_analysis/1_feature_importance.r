# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                Visualize Feature Importance
#
#                  Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(foreach)
library(doParallel)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(2)
ncore <- 4
registerDoParallel(cores = ncore)

barplot_r2_increase <- function(r2_df, age_include = TRUE) {
    # Map variable short names to full names
    variable_names <- c(
        age = "Age",
        sur_cover = "Surrounding Mature Forest Cover",
        num_fires = "Number of Fires",
        dist = "Distance to Nearest Mature Forest",
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
        mean_pdsi = "Mean Palmer Drought Severity Index",
        topography = "Topography",
        lu_sum_10 = "Pasture Years",
        lu_sum_20 = "Perennial Crop Years",
        lu_sum_30 = "Annual Crop Years",
        last_lu = "Last Land Use"
    )


    r2_df <- r2_df[r2_df$mean_r2_diff > 0.001, ]

    r2_df <- r2_df[order(r2_df$mean_r2_diff, decreasing = FALSE), ]

    r2_df$par <- factor(r2_df$par, levels = r2_df$par)

    custom_colors <- c(
        "#665191", "#a05195",
        "#d45087", "#f95d6a", "#ff7c43", "#ffa600",
        "#ffc300", "#ffda6a"
    )

    if (age_include) {
        custom_colors <- c("#003f5c", custom_colors)
    } else {
        r2_df <- r2_df[!(r2_df$par == "age"), ]
    }

    custom_colors <- rev(custom_colors[1:nrow(r2_df)])

    # # Interpolate gradient spanning all colors, matching number of categories
    # n_cats <- length(levels(r2_df$par))
    # color_palette <- colorRampPalette(custom_colors)(n_cats)

    p <- ggplot(r2_df, aes(x = group, y = mean_r2_diff, fill = par)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_y_continuous(
            name = "R² Increase",
        ) +
        scale_fill_manual(values = custom_colors, labels = variable_names[levels(r2_df$par)], name = "Variable") +
        theme_minimal(base_size = 12) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_blank(),
            legend.position = "right",
            text = element_text(size = 12),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks.y = element_line(color = "black")
        ) +
        scale_y_continuous(breaks = seq(0, 0.4, by = 0.10))

    if (!age_include) {
        p <- p + theme(
            axis.text.y = element_blank(),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    }

    return(p)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------- Amazon barplot (Nearest Mature) ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


r2_df <- read.csv("./0_results/0_r2_nearest_mature.csv")
r2_df$group <- "nearest_mature"

p <- barplot_r2_increase(r2_df, age_include = TRUE)
asymptote <- "nearest_mature"
ggsave(
    paste0("./0_results/figures/model_performance_", asymptote, ".jpg"),
    plot = p,
    width = 8,
    height = 10,
    dpi = 300
)

p <- barplot_r2_increase(r2_df, age_include = FALSE)
ggsave(
    paste0("./0_results/figures/model_performance_", asymptote, "_zoomed.jpg"),
    plot = p,
    width = 8,
    height = 10,
    dpi = 300
)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------- Atlantic Forest ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

r2_df <- read.csv("./0_results/0_r2_land_use_aggregated_all_1.csv")
r2_df$group <- "aggregated_all"

p <- barplot_r2_increase(r2_df, age_include = TRUE)

ggsave(
    paste0("./0_results/figures/model_performance_", "aggregated_all", ".jpg"),
    plot = p,
    width = 8,
    height = 10,
    dpi = 300
)


# ------------------------------------------------- #
# Figure - R2 per Asymptote with age_only
# ------------------------------------------------- #


r2_asymptote <- read.csv("./0_results/0_asymptotes.csv")
# Filter for the variables of interest
r2_asymptote <- r2_asymptote %>%
    filter(
        asymptote %in% c("nearest_mature", "quarter_biomass", "full_amazon"),
        basic_pars_name == "lag"
    )

p <- ggplot(
    r2_asymptote %>% mutate(asymptote = reorder(asymptote, mean_r2)),
    aes(x = asymptote, y = mean_r2)
) +
    geom_bar(stat = "identity", fill = "#003f5c") +
    geom_errorbar(
        aes(
            ymin = mean_r2 - sd_r2,
            ymax = mean_r2 + sd_r2
        ),
        width = 0.4,
        linewidth = 1,
        color = "black"
    ) +
    theme_minimal() +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    )

ggsave(
    paste0("./0_results/figures/figure_2_asymptote_barplot.jpg"),
    plot = p,
    width = 6,
    height = 10,
    dpi = 300
)
