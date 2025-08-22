# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#         Predicted vs Observed AGB from Satellite Data
#
#                    Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

library(tidyverse)
library(ggplot2)

set.seed(1)

data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)

# Fit the model on the full data
norm_data <- normalize_independently(data)
train_stats <- norm_data$train_stats
norm_data <- norm_data$train_data

pars_init <- find_combination_pars(
    basic_pars = basic_pars_options[["lag"]],
    data_pars = data_pars_options(colnames(data))[["all_mean_climate"]],
    norm_data
)

final_model <- run_optim(norm_data, pars_init, conditions)

pred <- growth_curve(final_model$par, data = norm_data, lag = final_model$par["lag"])

# Assuming pred = predicted AGB, obs = norm_data$biomass
df <- data.frame(
    Predicted = pred,
    Observed = norm_data$biomass
)

ext <- ggplot(df, aes(x = Predicted, y = Observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 2) +
    labs(
        x = "Predicted Biomass (Mg/ha)",
        y = "Observed Biomass (Mg/ha)"
    ) +
    coord_cartesian(expand = FALSE) +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black", size = 28, family = "Helvetica"),
        axis.text = element_text(color = "black", size = 18, family = "Helvetica"),
        legend.position = "none"
    )

# Save to file
ggsave("./0_results/figures/extended/predicted_vs_observed_satellite.png", 
    plot = ext
)
