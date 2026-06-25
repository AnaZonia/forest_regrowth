
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

data <- subset(data, ave(seq_along(ecoreg), ecoreg, FUN = length) >= 1000)


df_results <- data.frame()

for (ecoregion in unique(data$ecoreg)) {
    print(ecoregion)

    df_ecoreg <- subset(data, data$ecoreg == ecoregion)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(df_ecoreg))[["all"]]

    cv_results <- cross_validate(df_ecoreg, basic_pars, data_pars, conditions, folds = 3)

    result <- c(ecoregion, mean(cv_results[[3]]), sd(cv_results[[3]]))
    df_results <- rbind(df_results, result)
}

df_results

write_rds(df_results, file = "./0_results/lag_per_ecoregion.rds")
