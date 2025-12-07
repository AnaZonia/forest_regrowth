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


data <- import_data(paste0('grid_10k_amazon_secondary_modis'), biome = 1, n_samples = 240000, asymptote = "nearest_mature")

# remove mean_aet and rename modis_aet to mean_aet
data <- data %>%
    select(-mean_aet) %>%
    rename(mean_aet = modis_aet)

# age_one <- data %>% filter(age == 1)
# mean(age_one$biomass)

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[["all"]]

cv_all <- cross_validate(data, basic_pars, data_pars, conditions, 8)





lag_mean <- data.frame(
            mean_lag = mean(cv_all[[3]]),
            sd_lag = sd(cv_all[[3]])
        )
lag_mean

write.csv(lag_mean, file = "./0_results/0_lag.csv", row.names = FALSE)





write.csv(pars, file = "./0_results/0_parameters.csv", row.names = FALSE)

col_means <- colMeans(pars, na.rm = TRUE) # numeric named vector
means_df <- as.data.frame(t(col_means)) # one row, columns = predictors

# 2. Palette: red (min) -> white (0) -> blue (max)
pal <- col_numeric(
    palette = colorRampPalette(c("red", "white", "blue"))(100),
    domain  = range(col_means, na.rm = TRUE)
)

# 3. Assign a color per cell
cell_colors <- pal(as.numeric(means_df[1, ]))

datatable(round(means_df, 4), options = list(dom = "t")) %>%
    formatStyle(
        columns = colnames(means_df),
        backgroundColor = styleEqual(
            levels = as.numeric(means_df[1, ]),
            values = cell_colors
        )
    )










# variables to check
vars <- c("mean_aet", "mean_temp", "mean_pr", "mean_vpd", "phh2o", "mean_srad", "mean_def")

# 1) Correlation matrix (pairwise complete obs)
cor_mat <- cor(data[, vars], use = "pairwise.complete.obs")
print(cor_mat)

# Optional quick visualization
pairs(data[, vars])

# 2) VIFs from a simple linear model
#    (use any numeric response; here assume 'regrowth' exists in `data`)
library(car)

m <- lm(biomass ~ mean_temp + mean_vpd + phh2o + mean_srad + mean_pr + mean_aet,
    data = data
)

summary(m)

# summary(lm(biomass ~ mean_vpd, data = data))

vif_values <- car::vif(m)
print(vif_values)


colnames(data)

summary(lm(biomass ~ modis_aet, data = data))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------- R2 increase per Asymptote ------------------ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


groups <- c("full_amazon", "nearest_mature")
group_names <- c("Amazon-wide Average", "Nearest Neighbor")

results <- data.frame()

for (i in seq_along(groups)) {
    asymptote <- groups[i]
    data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 240000, asymptote = asymptote)

    basic_pars <- basic_pars_options[["lag"]]
    # data_pars <- data_pars_options(colnames(data))[["all"]]

    data_pars <- c("mean_aet", "sur_cover", "num_fires", "dist", "sand")

    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    cv_results

    write.csv(cv_results[[2]], file = paste0("./0_results/0_r2_", asymptote, "2.csv"), row.names = FALSE)

    result <- data.frame(
        asymptote = asymptote,
        mean_r2 = mean(cv_results[[1]]),
        sd_r2 = sd(cv_results[[1]])
    )
    results <- rbind(results, result)
    write.csv(results, file = "./0_results/0_asymptotes_all_pars2.csv", row.names = FALSE)
}

lm(biomass ~ mean_temp, data = data) %>% summary()


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

    cv_results <- cross_validate(data, basic_pars, data_pars, conditions)

    write.csv(cv_results[[2]], file = paste0("./0_results/0_r2_", land_use_aggregation, "_", biome, ".csv"), row.names = FALSE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------- Atlantic Forest ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


land_use_list <- list.files(
    path = "./0_data",
    pattern = "land_use",
    full.names = FALSE
)

# Step 2: Initialize results dataframe
results <- data.frame()

# Step 3: Aggregate results from each group and filter
for (land_use_aggregation in land_use_list) {
    # Read CSV for each aggregation group
    r2_df <- read.csv(
        paste0("./0_results/0_r2_", land_use_aggregation, "_4.csv")
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

# Step 4: Reshape to wide format, each group gets mean and sd columns
results_wide <- results %>%
    select(par, mean_r2_diff, sd_r2_diff, group) %>%
    pivot_wider(
        names_from = group,
        values_from = c(mean_r2_diff, sd_r2_diff)
    )

# Step 5: Set row names and remove parameter column
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

# Step 8: Print an xtable LaTeX table (include booktabs for quality lines)
print(
    xtable(table_combined),
    include.rownames = FALSE,
    booktabs = TRUE
)



