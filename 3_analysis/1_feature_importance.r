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
        topography = "Topography"
    )


    r2_df <- r2_df[r2_df$mean_r2_diff > 0.001, ]

    r2_df <- r2_df[order(r2_df$mean_r2_diff, decreasing = FALSE), ]

    r2_df$par <- factor(r2_df$par, levels = r2_df$par)

    custom_colors <- c("#665191", "#a05195",
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
        ) + scale_y_continuous(breaks = seq(0, 0.4, by = 0.10))

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------- Amazon barplot (two asymptotes) ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

groups <- c("full_amazon", "nearest_mature")
group_names <- c("Amazon-wide Average", "Nearest Neighbor")

results <- data.frame()

for (i in seq_along(groups)) {
    asymptote <- groups[i]
    data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 20000, asymptote = asymptote)
    norm_data <- normalize_independently(data)$train_data
    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all"]]
    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    write.csv(cv_results[[2]], file = paste0("./0_results/0_r2_", asymptote, ".csv"), row.names = FALSE)

    r2_df <- cv_results[[2]]
    r2_df$group <- asymptote
    p <- barplot_r2_increase(r2_df)

    result <- data.frame(
        asymptote = asymptote,
        mean_r2 = mean(cv_results[[1]]),
        sd_r2 = sd(cv_results[[1]])
    )
    results <- rbind(results, result)
    write.csv(results, file = "./0_results/0_asymptotes_all_pars.csv", row.names = FALSE)

    ggsave(
        paste0("./0_results/figures/model_performance_", asymptote, ".jpg"),
        plot = p,
        width = 8,
        height = 10,
        dpi = 300
    )
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------- Atlantic Forest ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



land_use_list <- list.files(paste0("./0_data"), pattern = "land_use", full.names = FALSE)

biomes_r2 <- list()
biomes_par_importance <- list()

results <- data.frame()

for (land_use_aggregation in land_use_list) {
    all_results <- list()
    group_names <- c("Amazon", "Atlantic Forest")

    for (biome in c(1, 4)) {
        data <- import_data(land_use_aggregation, biome = biome, n_samples = 10000, asymptote = "nearest_mature")
        basic_pars <- basic_pars_options[["lag"]]
        data_pars <- data_pars_options(colnames(data))[["all"]]

        cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

        write.csv(cv_results[[2]], file = paste0("./0_results/0_r2_", land_use_aggregation, "_", biome, ".csv"), row.names = FALSE)


        result <- data.frame(
            biome = biome,
            land_use_aggregation = land_use_aggregation,
            mean_r2 = mean(cv_results[[1]]),
            sd_r2 = sd(cv_results[[1]])
        )
        results <- rbind(results, result)
        write.csv(results, file = "./0_results/0_land_use_R2.csv", row.names = FALSE)

    }

}

# ------------------------------------------------- #
# Figure - R2 per Asymptote with age_only
# ------------------------------------------------- #


R2_asymptote <- read.csv("./0_results/0_asymptotes.csv")
# Filter for the variables of interest
R2_asymptote <- R2_asymptote %>%
    filter(asymptote %in% c("nearest_mature", "quarter_biomass", "full_amazon"))

p <- ggplot(
    R2_asymptote %>% mutate(asymptote = reorder(asymptote, mean_r2)),
    aes(x = asymptote, y = mean_r2)
) +
    geom_bar(stat = "identity", fill = "#003f5c") +
    coord_flip() +
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
    width = 10,
    height = 6,
    dpi = 300
)