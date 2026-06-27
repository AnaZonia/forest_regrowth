
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
source("2_modelling/3_genetic_algorithm.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- One lag per ecoregion of the Amazon ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# choose which ecoregions to include based on the number of secondary forest pixels 

data <- import_data("ecoreg_stratify", biome = 1, n_samples = 240000, asymptote = "nearest_mature", categorical = c("topography", "last_lu"))

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

norm_data <- normalize_independently(data)
train_data <- norm_data$train_data
train_stats <- norm_data$train_stats

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[["all"]]

pars_init <- forward_selection(basic_pars, data_pars, train_data)


df_results <- data.frame()

for (ecoregion in unique(data$ecoreg)) {
    print(ecoregion)

    df_ecoreg <- subset(data, data$ecoreg == ecoregion)

    df_ecoreg <- apply_min_max_scaling(df_ecoreg, train_stats)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(df_ecoreg))[["all"]]

    min_pars <- optim_ga(ini_par, df_ecoreg)

    df_results <- rbind(df_results, cbind(min_pars, ecoregion))
    print(df_results)
    write_rds(df_results, file = "./0_results/lag_per_ecoregion.rds")
}

# ------


# data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 30000, asymptote = "nearest_mature", categorical = categorical)

# # identify ecoregion columns
# eco_cols <- grep("^ecoreg_", names(data), value = TRUE)

# # count rows with value 1 in each ecoregion column
# eco_counts <- colSums(data[eco_cols] == 1, na.rm = TRUE)
# eco_counts

# # keep only ecoregion columns with at least 1000 rows equal to 1
# keep_eco <- names(eco_counts)[eco_counts >= 1000]

# # remove the rest
# data2 <- data[, c(setdiff(names(data), eco_cols), keep_eco)]


# data2_481 <- subset(data2, data2$ecoreg_481 == 1)

# norm_data <- normalize_independently(data2)
# train_data <- norm_data$train_data
# train_stats <- norm_data$train_stats

# basic_pars <- basic_pars_options[["lag"]]
# data_pars <- data_pars_options(colnames(train_data))[["all"]]

# pars_init <- forward_selection(basic_pars, data_pars, train_data)

# pars_init

# # 32.5

# ini_par <- pars_init[[1]]

# for (j in 2:ncore) {
#     for (name in names(ini_par)) {
#         ini_par[j, name] <- c(
#             ini_par[1, name] * (1.5 * runif(1) + .5)
#         )
#     }
# }

# ini_par

# par = ini_par

# norm_data = train_data

# min_pars <- optim_ga(ini_par, train_data)





# what if I keep everything the same, changing only the lag?
# use all of the same fit parameters used for the full Amazon, but fit the lag to each specific ecoregion?