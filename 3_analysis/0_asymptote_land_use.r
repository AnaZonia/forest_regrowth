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



print(head(iter_df))

data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 50000, asymptote = "nearest_mature")

data_pars_name <- "all"

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

pars_init <- find_combination_pars(basic_pars, data_pars, data)



# Identify ecoreg rows
ecoreg_mask <- grepl("^ecoreg_", df$par)

# Check if all ecoreg r2 values are equal
if (length(unique(df$r2[ecoreg_mask])) == 1) {
    # Get the unique r2 value
    ecoreg_r2 <- unique(df$r2[ecoreg_mask])
    # Remove all ecoreg_* rows
    df <- df[!ecoreg_mask, ]
    # Add single "ecoreg" row
    df <- rbind(df, data.frame(par = "ecoreg", r2 = ecoreg_r2))
    # Optional: reorder rows if needed
    rownames(df) <- NULL
}

print(df)

cv_results <- cross_validate(data, basic_pars, data_pars, conditions)
cv_results








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
        data <- import_data("grid_10k_amazon_secondary", biome_num = 1, n_samples = 50000, asymptote = asymptote)

        data_pars_name <- "age_only"

        basic_pars <- basic_pars_options[[basic_pars_name]]
        data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

        cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

        result <- data.frame(
            basic_pars_name = basic_pars_name,
            asymptote = asymptote,
            mean_r2 = mean(cv_results),
            sd_r2 = sd(cv_results)
        )

        print(result)
        results <- rbind(results, result)
        write.csv(results, file = "./0_results/0_asymptotes.csv", row.names = FALSE)
    }
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

results <- data.frame()

land_use_list <- list.files(paste0("./0_data"), pattern = "land_use", full.names = FALSE)

for (biome in c(1, 4)) {
    print(biome)

    for (land_use_aggregation in land_use_list) {

        data <- import_data(land_use_aggregation, biome_num = biome, n_samples = 50000)

        for (data_pars_name in c("age_only", "land_use", "fires")) {

            basic_pars <- basic_pars_options[["lag"]]
            data_pars <- data_pars_options(colnames(data))[[data_pars_name]]

            cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

            result <- data.frame(
                land_use_aggregation = land_use_aggregation,
                data_pars_name = data_pars_name,
                biome = biome,
                mean_r2 = mean(cv_results),
                sd_r2 = sd(cv_results)
            )

            print(result)
            results <- rbind(results, result)
            write.csv(results, file = "./0_results/0_land_use.csv", row.names = FALSE)
        }
    }
}