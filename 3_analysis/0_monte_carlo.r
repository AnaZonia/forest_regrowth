
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

# run monte carlo with sd from biomass as distribution of error (error propagation)

data <- import_data("monte_carlo", biome = 1, n_samples = 30000, asymptote = "nearest_mature")

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[["all"]]

error_prop_results <- error_prop(data, basic_pars, data_pars, conditions)

error_prop_results

# get the distribution of lag values


cv_results <- cross_validate(data, basic_pars, data_pars, conditions, 5)

mean(cv_results[[1]])

mean(cv_results[[3]])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- One lag per ecoregion of the Amazon ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

csv_files <- list.files(paste0("./0_data/monte_carlo"), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()

df <- subset(df, ave(seq_along(ecoreg), ecoreg, FUN = length) >= 500)

for (ecoregion in unique(df$ecoreg)) {
    df_ecoreg <- subset(df, df$ecoreg == ecoregion)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all"]]

    cv_results <- cross_validate(data, basic_pars, data_pars, conditions, 5)

    result <- data.frame(
        mean_r2 = mean(cv_results[[1]]),
        mean_lag = mean(cv_results[[3]])
    )
    print(result)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- Future predictions - Bezerra map ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# make future projections with the ecoregion-specific models or with the general model from error propagation?

