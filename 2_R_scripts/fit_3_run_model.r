# Compares results for the fits with:

# 2. Different land use aggregations

library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR) # For PCA
library(factoextra) # Optional, for PCA visualization
library(cluster) # For k-means clustering

# Source external R scripts for data import and function definitions
source("./2_R_scripts/fit_1_import_data.r")
source("./2_R_scripts/fit_2_functions.r")

# Set up parallel processing
set.seed(1)
ncores <- 35
registerDoParallel(cores = ncores)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# List of climatic parameters that change yearly
climatic_pars <- c("srad", "soil", "temp", "vpd", "aet", "def", "pdsi", "pr")
land_use <- c("lu", "fallow", "num_fires")
landscape <- c("dist", "sur_cover", "nearest_mature_biomass")
categorical <- c("ecoreg", "topography") # , "last_lu")
# indig, protec and floodable_forests are already boolean

# Define conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0', 'pars["k0"] < 0')

non_data_pars <- c("k0", "B0", "theta", "lag")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Switches ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# basic_pars_name <- "lag"
basic_pars_name <- "intercept"
# basic_pars_name <- "intercept_age_multiplicative"

# interval <- "5yr", "10yr", "15yr", "all"
# import_name <- "aggregated"
# export_name <- import_name

biome <- 1

# data_pars_names <- "land_use_landscape_only"
data_pars_names <- "all"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dataframe <- import_data("./0_data/unified_fc.csv", biome = biome, n_samples = 15000)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define basic parameter sets for modeling
if (basic_pars_name == "lag") {
    basic_pars <- c("lag", "k0", "theta")
} else if (basic_pars_name == "intercept") {
    basic_pars <- c("k0", "B0", "theta")
    # age, unnormalized, multiplies over the entire linear combination k
} else if (basic_pars_name == "intercept_age_multiplicative") {
    basic_pars <- c("k0", "B0", "theta", "age")
    # age is just one more parameter, treated like the others
    # just for sanity check, should be pretty much the same as intercept
}

colnames <- colnames(dataframe)

exclusion_pattern <- paste(c("age", "biomass", "nearest_biomass", paste0(climatic_pars, "_")), collapse = "|")

if (data_pars_names == "land_use_landscape_only") {
    data_pars = c(colnames[grepl(paste0(c(land_use, landscape), collapse = "|"), colnames)])
} else if (data_pars_names == "all") {
    data_pars = colnames[!grepl(paste0(c(exclusion_pattern, categorical, landscape, land_use), collapse = "|"), colnames)]
}

c(colnames_filtered[!grepl(paste0(c(land_use, landscape, "protec", "indig", "nitro", "sand"), collapse = "|"), colnames_filtered)]), #only climatic
c(colnames_filtered[!grepl(paste0(c(land_use, landscape), collapse = "|"), colnames_filtered)]), #climatic, ecoregion, topography, soil type
c(colnames_filtered[!grepl(paste0(land_use, collapse = "|"), colnames_filtered)]), #climatic, distance, sur_cover
c(colnames_filtered[!grepl(paste0(land_use, collapse = "|"), colnames_filtered)], 'num_fires'), #climatic, distance, sur_cover
c(colnames_filtered[!grepl("ecoreg", colnames_filtered)])
# c(colnames_filtered_no_mean_climate, climatic_pars)



print(data_pars)
print(basic_pars)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----------------------------- K-Fold Cross-Validation ---------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Function Description:
#   This function performs k-fold cross-validation (with k=5) on the provided data.
#   It splits the data into five folds, uses each fold as a test set once while training
#   on the remaining folds, and then applies the specified `run_function` to train the model
#   and calculate the R-squared values for each fold.
#
# Arguments:
#   data        : The full dataset to be split into training and test sets for cross-validation.
#   run_function: The function used to train the model and generate predictions (either "run_optim", "run_lm")
#   pars_iter   : The parameters to be passed to `run_function` for model training.
#   conditions  : Additional conditions or constraints to be passed to `run_function`
#                 (optional, default is NULL).
#
# Returns:
#   A list with the following elements:
#     - `r2`     : The mean R-squared value across all folds.
#     - `r2_sd`  : The standard deviation of the R-squared values across all folds.
#     - `pars`    : The model parameters corresponding to the fold with the highest R-squared value.
#
# Notes:
#   - The function uses random sampling to assign data points to each of the five folds.
#   - The function assumes that `run_function` takes as input the training data, parameters,
#     conditions, and test data, and returns the model output.
#   - The R-squared values from each fold are stored in `r2_list`, and the best model
#


indices <- sample(c(1:5), nrow(dataframe), replace = TRUE)
dataframe$pred_cv <- NA
dataframe$pred_final <- NA
r2_list <- numeric(5)


for (index in 1:5) {
    # Define the test and train sets
    test_data <- dataframe[indices == index, -grep("pred", names(dataframe))]
    train_data <- dataframe[indices != index, -grep("pred", names(dataframe))]
    # Normalize training and test sets independently, but using training data's min/max for both
    norm_data <- normalize_independently(train_data, test_data)

    train_data <- norm_data$train_data
    test_data <- norm_data$test_data

    # Function to perform direct optimization
    pars_init <- find_combination_pars(basic_pars, data_pars, train_data)

    # Run the model function on the training set and evaluate on the test set
    model <- run_optim(train_data, pars_init, conditions)
    pred_cv <- growth_curve(model$par, test_data)

    # save the predicted values of each iteration of the cross validation.
    dataframe$pred_cv[indices == index] <- pred_cv
    r2 <- calc_r2(dataframe[indices == index, ], pred_cv)
    r2_list[index] <- r2
    print(r2)
}


print(mean(r2_list))













# # Function to perform clustered optimization
# # can cluster by anthropogenic, environmental, or both.
# clustered_optimization <- function(iterations, n_clusters = 5, cluster_by = "anthro") {
#     for (iter in 1:nrow(iterations)) {
#         r2_cluster <- numeric(n_clusters)

#         i <- iterations$interval[iter]
#         j <- iterations$data_par[iter]
#         k <- iterations$biome[iter]
#         l <- iterations$basic_par[iter]

#         data_iter <- dataframes[[i]][[k]]
#         data_pars_name <- data_pars_names[[j]]
#         biome_name <- biomes[[k]]
#         basic_pars_name <- basic_pars_names[[l]]

#         data_iter$pred <- NA
#         if (cluster_by == "enviro") {
#             columns_cluster <- names(data_iter)[!grepl(paste0(c(land_use, landscape), collapse = "|"), names(data_iter))]
#         } else if (cluster_by == "anthro") {
#             columns_cluster <- names(data_iter)[grepl(paste0(c(land_use, landscape), collapse = "|"), names(data_iter))]
#         }

#         data_norm <- data_iter[, names(data_iter) %in% columns_cluster]
#         data_norm <- normalize_independently(data_norm)$train_data

#         pca <- PCA(data_norm, ncp = n_clusters, graph = FALSE)
#         data_norm_pca <- as.data.frame(pca$ind$coord)
#         kmeans_result <- kmeans(data_norm_pca, centers = n_clusters, nstart = 20)
#         data_iter$cluster <- kmeans_result$cluster

#         for (cluster_id in 1:n_clusters) {
#             data_cluster <- subset(data_iter, cluster == cluster_id)
#             pars_iter <- find_combination_pars(iter, data_cluster)
#             optim_cv_output <- cross_valid(data_cluster, pars_iter, conditions)
#             r2_cluster[cluster_id] <- optim_cv_output$r2_final
#         }
#         # c(mean(r2_cluster), sd(r2_cluster))
#     }

#     return(c(mean(r2_cluster), sd(r2_cluster)))
# }


# if (run_mode == "direct") {
#     direct_optimization(iterations)
# } else if (run_mode == "clustered") {
#     clustered_optimization(iterations, n_clusters = 5, cluster_by = cluster_by)
# }


# Export lat lon 
