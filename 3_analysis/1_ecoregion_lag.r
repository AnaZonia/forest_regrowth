
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#      Compare R2 of different asymptotes and land use aggregations
#
#                     Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(foreach)
library(doParallel)
library(tidyverse)
library(xtable)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- One lag per ecoregion of the Amazon ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 10k per ecoregion
# check range and distribution per ecoregions

data <- import_data("ecoreg_stratify", biome = 1, n_samples = 240000, asymptote = "nearest_mature")

results <- c()

# error_prop <- readRDS("./0_results/0_error_prop.rds")

apply_min_max_scaling <- function(data, train_stats) {
    # Apply Min-Max scaling to each variable in the data
    for (i in seq_along(train_stats$variable)) {
        var <- train_stats$variable[i]
        print(var)
        data[[var]] <- (data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])
    }
    return(data)
}


for (ecoregion in c(508, 518, 507, 481, 476, 497)) {
    print(ecoregion)

    df_ecoreg <- subset(data, data$ecoreg == ecoregion)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(df_ecoreg))[["all"]]

    cv_results <- cross_validate(df_ecoreg, basic_pars, data_pars, conditions, folds = 2)

    # norm_df <- apply_min_max_scaling(future[, names(future) %in% train_stats$variable], train_stats)
    # print(cv_results[[3]])

    results <- c(results, mean(cv_results[[3]]))
}

nrow_ecoreg <- c()
for (ecoregion in unique(df$ecoreg)) {
    df_ecoreg <- subset(df, df$ecoreg == ecoregion)
    nrow_ecoreg <- c(nrow_ecoreg, nrow(df_ecoreg))
}

df_results <- data.frame(ecoreg = unique(df$ecoreg), lag = results, nrow = nrow_ecoreg)

df_results

