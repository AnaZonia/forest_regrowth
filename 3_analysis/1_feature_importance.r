# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                Visualize Feature Importance
#
#                  Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(ggplot2)
library(tidyverse)
library(RColorBrewer)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")
source("2_modelling/2_permutation_importance.r")

# Set up parallel processing
set.seed(2)
ncore <- 4
registerDoParallel(cores = ncore)

barplot_r2_increase <- function(all_results, group_names) {

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
        mean_pdsi = "Mean Palmer Drought Severity Index"
    )

    r2_df$asymptote <- "nearest_mature"

    r2_df$par <- factor(r2_df$par, levels = r2_df$par)

    r2_df <- r2_df %>%
        filter(par != "age")

    custom_colors <- c(
        "#003f5c", "#2f4b7c", "#665191", "#a05195",
        "#d45087", "#f95d6a", "#ff7c43", "#ffa600", "#ffc300", "#ffda6a"
    )

    recycled_colors <- rep(custom_colors, length.out = 10)

    p <- ggplot(r2_df, aes(x = asymptote, y = r2_diff, fill = par)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_y_continuous(
            name = "R² Increase",
        ) +
        scale_fill_manual(values = recycled_colors, labels = variable_names[levels(r2_df$par)], name = "Variable") +
        theme_minimal(base_size = 12) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_blank(),
            legend.position = "right",
            text = element_text(size = 12),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16, face = "bold"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()
        )

    p


    return(p)
}
ggplot(r2_df, aes(x = par, y = r2_diff, fill = group)) +
    geom_bar(stat = "identity") +
    labs(y = "R² Increase", x = "Variable") +
    theme_minimal()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------- Amazon barplot (two asymptotes) ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

groups <- c("full_amazon", "nearest_mature")
group_names <- c("Amazon-wide Average", "Nearest Neighbor")

for (i in seq_along(groups)) {
    asymptote <- groups[i]
    data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000, asymptote = asymptote)
    norm_data <- normalize_independently(data)$train_data
    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
    r2 <- calc_r2(norm_data, pred)
    importance_results <- calculate_permutation_importance(model, norm_data, data_pars)
    importance_results$importance_scaled <- importance_results$importance_pct * r2 / 100
    importance_results$group <- group_names[i]


}

p <- barplot_feat_importance(all_results, group_names)
p
ggsave(
    "./0_results/figures/model_performance_def_phh2o_2.jpg",
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------- Atlantic Forest ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



land_use_list <- list.files(paste0("./0_data"), pattern = "land_use", full.names = FALSE)

for (land_use_aggregation in land_use_list) {
    all_results <- list()
    group_names <- c("Amazon", "Atlantic Forest")

    i <- 1
    for (biome in c(1, 4)) {
        data <- import_data(land_use_aggregation, biome_num = biome, n_samples = 10000, asymptote = "nearest_mature")
        norm_data <- normalize_independently(data)$train_data
        basic_pars <- basic_pars_options[["lag"]]
        data_pars <- data_pars_options(colnames(data))[["all_mean_climate"]]
        init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
        model <- run_optim(norm_data, init_pars, conditions)
        pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
        r2 <- calc_r2(norm_data, pred)
        importance_results <- calculate_permutation_importance(model, norm_data, data_pars)
        importance_results$importance_scaled <- importance_results$importance_pct * r2 / 100
        importance_results$group <- group_names[i]
        all_results[[i]] <- importance_results
        i <- i + 1
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