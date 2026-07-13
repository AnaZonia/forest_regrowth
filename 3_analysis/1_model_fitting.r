
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- Fit model through error propagation ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   - Obtain mean lag and R² through uncertainty propagation
#     - error_prop.rds
#     - r2_full_amazon_error_prop.csv
#   - Obtain R² through different asymptote aggregations
#     - asymptotes.csv
#   - Compare R² increase from each Land Use aggregation
#   - Compare R² with land use as predictors for the Atlantic Forest
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(foreach)
library(doParallel)
library(tidyverse)
library(xtable)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_error_propagation.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore = 4
registerDoParallel(cores = ncore)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------- Lag and R² - uncertainty propagation ---------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("grid_10k_amazon_uncertainty_propagation", biome = 1, n_samples = 30000, asymptote = "nearest_mature", categorical = categorical)


basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[["all"]]

error_prop_results <- error_prop(data, basic_pars, data_pars, conditions, 1000)

results <- data.frame(
    mean_r2 = mean(error_prop_results[[1]]),
    sd_r2 = sd(error_prop_results[[1]]),
    mean_lag = mean(error_prop_results$pars[["lag"]]),
    sd_lag = sd(error_prop_results$pars[["lag"]])
)

results

# write_rds(error_prop_results, file = "./0_results/error_prop.rds")
# write_csv(results, "./0_results/r2_full_amazon_error_prop.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- Asymptote Comparisons ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Compare model R² performance for each asymptote reference:
#   - nearest_mature
#   - ecoreg_biomass
#   - quarter_biomass
#   - full_amazon
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

results <- data.frame()

for (asymptote in c("nearest_mature", "quarter_biomass", "full_amazon")) {
    data <- import_data("grid_10k_amazon_uncertainty_propagation", biome = 1, n_samples = 30000, asymptote = asymptote, categorical = categorical)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["age_only"]]

    error_prop_results <- error_prop(data, basic_pars, data_pars, conditions, 500)

    result <- data.frame(
        asymptote = asymptote,
        mean_r2 = mean(error_prop_results[[1]]),
        sd_r2 = sd(error_prop_results[[1]])
    )

    print(result)
    results <- rbind(results, result)
    write.csv(results, file = "./0_results/asymptotes.csv", row.names = FALSE)
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
    biome = 4

    data <- import_data(land_use_aggregation, biome = biome, n_samples = 30000, asymptote = "nearest_mature")

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(data))[["all"]]

    cv_results <- error_prop(data, basic_pars, data_pars, conditions, 1000)

    write.csv(cv_results[[2]], file = paste0("./0_results/r2_", land_use_aggregation, "_", biome, ".csv"), row.names = FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------- Atlantic Forest ------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


land_use_list <- list.files(
    path = "./0_data",
    pattern = "land_use",
    full.names = FALSE
)

results <- data.frame()

# Aggregate results from each group and filter
for (land_use_aggregation in land_use_list) {
    # Read CSV for each aggregation group
    r2_df <- read.csv(
        paste0("./0_results/r2_", land_use_aggregation, "_4.csv")
    )
    r2_df$group <- land_use_aggregation # Tag with group

    # Filter for relevant r2 differences
    r2_df <- r2_df[r2_df$mean_r2_diff > 0.001, ]
    # Optional: Order results by descending mean difference
    r2_df <- r2_df[order(r2_df$mean_r2_diff, decreasing = TRUE), ]

    # Append to cumulative results
    results <- rbind(results, r2_df)
}

# Keep naming consistent for aggregated and non_aggregated
results$par[results$par == "lu_sum_15"] <- "lu_sum_10"

# Reshape to wide format, each group gets mean and sd columns
results_wide <- results %>%
    select(par, mean_r2_diff, sd_r2_diff, group) %>%
    pivot_wider(
        names_from = group,
        values_from = c(mean_r2_diff, sd_r2_diff)
    )

# Set row names and remove parameter column
results_wide <- as.data.frame(results_wide)
rownames(results_wide) <- results_wide$par
results_wide$par <- NULL


# Remove "land_use_" prefix for clearer labels
col_labels <- sub("land_use_", "", land_use_list)

# Dynamically build the column names to select ("mean \\pm sd" columns for each group)
selected_cols <- sapply(col_labels, function(lab) lab)

# Generate combined columns (if not already made), using your previous loop
for (i in seq_along(land_use_list)) {
    mean_col <- paste0("mean_r2_diff_", land_use_list[i])
    sd_col <- paste0("sd_r2_diff_", land_use_list[i])
    out_col <- col_labels[i]
    results_wide[[out_col]] <- sprintf("%.5f ± %.5f", results_wide[[mean_col]], results_wide[[sd_col]])
}

# Assemble the table with parameter descriptions and formatted columns
table_combined <- results_wide[, c(selected_cols)]
table_combined$par <- rownames(results_wide) # assumes par names are already readable, otherwise remap here
table_combined <- table_combined[, c("par", selected_cols)]

# Replace NA/missing with "---"
table_combined[is.na(table_combined)] <- "---"

# Rename parameters using a named vector
par_labels <- c(
    sur_cover = "Surrounding Mature Forest Cover",
    mean_aet = "Mean Actual Evapotranspiration",
    age = "Stand Age",
    num_fires = "Number of Fires",
    dist = "Distance to Nearest Mature Forest",
    topography = "Topography",
    lu_sum_10 = "Number of Pasture Years",
    phh2o = "Soil pH"
)

# Replace par codes with their full descriptions
table_combined$par <- par_labels[table_combined$par]
# Set empty string for parameter column name and readable labels for others
colnames(table_combined) <- c("", "Aggregated All", "Non-Aggregated 10yr", "Non-Aggregated 5yr", "Non-Aggregated All")
table_combined[is.na(table_combined) | table_combined == "NA ± NA"] <- "---"

print(
    xtable(table_combined),
    include.rownames = FALSE,
    booktabs = TRUE
)




