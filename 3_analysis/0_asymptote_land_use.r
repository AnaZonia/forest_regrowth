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

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore = 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Asymptote Comparisons ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Compare model R² performance for each asymptote reference:
#   - nearest_mature
#   - ecoreg_biomass
#   - quarter_biomass
#   - full_amazon)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


results <- data.frame()

for (asymptote in c("nearest_mature", "ecoreg_biomass", "quarter_biomass", "full_amazon")) {
    for (basic_pars_name in c("intercept", "lag")) {

        data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 30000, asymptote = asymptote)

        data_pars_name <- "age_only"

        basic_pars <- basic_pars_options[[basic_pars_name]]
        data_pars <- data_pars_options(colnames(data))[[data_pars_name]]    

        cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

        result <- data.frame(
            basic_pars_name = basic_pars_name,
            asymptote = asymptote,
            mean_r2 = mean(cv_results[[1]]),
            sd_r2 = sd(cv_results[[1]])
        )

        print(result)
        results <- rbind(results, result)
        write.csv(results, file = "./0_results/0_asymptotes.csv", row.names = FALSE)
    }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Average Lag expected ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 150000, asymptote = "nearest_mature")

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[["all"]]

cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

result <- data.frame(
    mean_lag = mean(cv_results[[3]]),
    sd_lag = sd(cv_results[[3]])
)

write.csv(result, file = "./0_results/0_lag.csv", row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- R2 increase per Asymptote ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


groups <- c("full_amazon", "nearest_mature")
group_names <- c("Amazon-wide Average", "Nearest Neighbor")

results <- data.frame()

for (i in seq_along(groups)) {
    asymptote <- groups[i]
    data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 30000, asymptote = asymptote)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all"]]
    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    write.csv(cv_results[[2]], file = paste0("./0_results/0_r2_", asymptote, ".csv"), row.names = FALSE)

    result <- data.frame(
        asymptote = asymptote,
        mean_r2 = mean(cv_results[[1]]),
        sd_r2 = sd(cv_results[[1]])
    )
    results <- rbind(results, result)
    write.csv(results, file = "./0_results/0_asymptotes_all_pars.csv", row.names = FALSE)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Land Use Comparisons ------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Compare R² increase from each Land Use aggregation:
#   - non_aggregated_all
#   - aggregated_all
#   - non_aggregated_5yr
#   - non_aggregated_10yr)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

land_use_list <- list.files(paste0("./0_data"), pattern = "land_use", full.names = FALSE)

biomes_r2 <- list()
biomes_par_importance <- list()

results <- data.frame()

for (land_use_aggregation in land_use_list) {
    all_results <- list()
    group_names <- c("Atlantic Forest") # "Amazon",

    for (biome in c(4)) { # 1
        data <- import_data(land_use_aggregation, biome = biome, n_samples = 30000, asymptote = "nearest_mature")

        basic_pars <- basic_pars_options[["lag"]]
        data_pars <- data_pars_options(colnames(data))[["all"]]

        cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

        write.csv(cv_results[[2]], file = paste0("./0_results/0_r2_", land_use_aggregation, "_", biome, ".csv"), row.names = FALSE)
    }
}

land_use_aggregation <- land_use_list[3]
