# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#         Predicted vs Observed AGB from Satellite Data
#
#                    Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

library(tidyverse)

set.seed(1)

data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 10000)

# remove the columns in data of the form ecoreg_xxx that have the value xxx not present in field_non_repeats$ecoreg
remove_unused_onehot <- function(data, field_non_repeats, prefix) {
    cols <- grep(paste0("^", prefix, "_"), colnames(data), value = TRUE)
    valid <- paste0(prefix, "_", unique(field_non_repeats[[prefix]]))
    cols_to_remove <- setdiff(cols, valid)
    data %>% select(-all_of(cols_to_remove))
}

data <- data %>%
    remove_unused_onehot(field_non_repeats, "ecoreg") %>%
    remove_unused_onehot(field_non_repeats, "topography")


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

# save plot as png
png("./0_results/figures/extended/predicted_vs_observed_satellite.png", width = 800, height = 600)
plot(pred, norm_data$biomass, xlab = "Predicted AGB", ylab = "Observed AGB", main = "Predicted vs Observed AGB")
abline(0, 1, col = "red", lty = 2)
dev.off()
