
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR) # For PCA
library(factoextra) # Optional, for PCA visualization
library(cluster) # For k-means clustering

# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

# Set up parallel processing
set.seed(1)
ncores <- 40
registerDoParallel(cores = ncores)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import_name <- "aggregated"
export_name <- import_name

# number of rows to be included in analysis
n_samples <- 10000

# List of climatic parameters that change yearly
climatic_pars <- c("srad", "soil", "temp", "vpd")
categorical <- c("ecoreg", "topography", "last_LU")
land_use <- c("lu", "fallow", "num_fires")
landscape <- c("distance", "sur_cover", "nearest_mature_biomass")

# Define conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0', 'pars["k0"] < 0')
non_data_pars <- c("k0", "B0", "theta", "lag")

# two biomes
# four dataframes


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define land-use history intervals to import four dataframes
# intervals <- list("5yr", "10yr", "15yr", "all")
intervals <- list("all")

# datafiles <- paste0("./new_data/", name_import, "_", intervals, ".csv")
# dataframes <- lapply(datafiles, import_data, convert_to_dummy = TRUE, process_climatic = FALSE)

dataframe <- read.csv("~/Documents/data/mapbiomas_heinrich_field.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Identify common columns across all dataframes of the same biome (across intervals)
# (avoids errors for parametesrs like last_LU, that not all dataframes may contain depending on how many rows were sampled)
# and filter out non-predictors
unique_colnames <- c()
data_pars <- c()

colnames_lists <- lapply(dataframes, function(df_list) {
    colnames(df_list[[i]])
})

# Step 2: Find the intersection of all column names
colnames_intersect <- Reduce(intersect, colnames_lists)
unique_colnames <- union(unique_colnames, colnames_intersect)
exclusion_pattern <- paste(c("age", "biomass", "nearest_mature_biomass", paste0(climatic_pars, "_")), collapse = "|")
# Filter colnames based on this pattern
colnames_filtered <- colnames_intersect[!grepl(exclusion_pattern, colnames_intersect)]
colnames_filtered_no_mean_climate <- colnames_filtered[!grepl(paste(paste0("mean_", climatic_pars), collapse = "|"), colnames_filtered)]

biome_pars <- list(
    c(colnames_filtered[grepl(paste0(c(land_use, landscape), collapse = "|"), colnames_filtered)])
)


# biome_pars <- list(
#     c(colnames_filtered[!grepl(paste0(c(land_use, landscape), collapse = "|"), colnames_filtered)]),
#     c(colnames_filtered[!grepl(paste0(land_use, collapse = "|"), colnames_filtered)]),
#     c(colnames_filtered[!grepl(paste0(landscape, collapse = "|"), colnames_filtered)]),
#     c(colnames_filtered),
#     c(colnames_filtered_no_mean_climate, climatic_pars)
# )

data_pars[[i]] <- biome_pars


data_pars_names <- c(
    "land_use"
)

# data_pars_names <- c(
#     "only_land_use_no_landscape",
#     "no_land_use",
#     "no_landscape",
#     "all",
#     "all_yearly_clim"
# )

# Define basic parameter sets for modeling
basic_pars <- list(
    # c("age", "k0", "B0"), #age is fit as its own independent predictor
    c("k0", "B0"), # k is multiplied by the age column
    c("lag", "k0")
)

basic_pars_names <- c(
    # "intercept_age",
    "intercept",
    "lag"
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

iterations <- expand.grid(
    basic_par = seq_along(basic_pars),
    data_par = seq_along(data_pars[[1]])
)

# basic_pars_fit_age <- which(sapply(basic_pars, function(x) "age" %in% x))
# data_pars_with_climatic <- which(sapply(data_pars[[1]], function(x) any(climatic_pars %in% x)))
# Remove rows where both conditions are met
# iterations <- iterations %>% filter(!(basic_par %in% basic_pars_fit_age & data_par %in% data_pars_with_climatic))
# iterations

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Function to perform direct optimization
direct_optimization <- function(iterations, batch_size = 20) {
    pars_iter_list <- list()
    for (iter in 1:nrow(iterations)) {
        i <- iterations$interval[iter]
        data <- dataframes[[i]]
        pars_iter <- find_combination_pars(iter, data)
        pars_iter_list[[iter]] <- pars_iter
        saveRDS(pars_iter_list, "./results/temp_pars_iter_list.rds")
    }

    num_batches <- ceiling(nrow(iterations) / batch_size)
    csv_files <- vector("list", num_batches)

    for (batch in 1:num_batches) {
        start <- (batch - 1) * batch_size + 1
        end <- min(batch * batch_size, nrow(iterations))

        results <- foreach(iter = start:end, .combine = rbind, .packages = "tidyverse") %dopar% {
            i <- iterations$interval[iter]
            j <- iterations$data_par[iter]
            k <- iterations$biome[iter]
            l <- iterations$basic_par[iter]

            pars_iter <- pars_iter_list[[iter]]
            data <- dataframes[[i]][[k]]
            pars_names <- data_pars_names[[j]]
            biome_name <- biomes[[k]]
            basic_pars_name <- basic_pars_names[[l]]

            cv_output <- cross_valid(data, pars_iter, conditions)
            row <- process_row(cv_output, basic_pars_name, intervals[[i]], pars_names, biome_name)

            if (!any(climatic_pars %in% names(pars_iter))) {
                data_lm <- dataframes_lm[[i]][[k]]
                cv_output <- cross_valid(data_lm, pars_iter)
                row_lm <- process_row(cv_output, "lm", intervals[[i]], pars_names, biome_name)
                row <- rbind(row, row_lm)
            }
            row
        }

        batch_file <- paste0("./new_data_yearly/", export_name, "_batch_", batch, "_results.csv")
        csv_files[[batch]] <- batch_file
        write.csv(results, batch_file, row.names = FALSE)
    }

    all_results <- do.call(rbind, lapply(csv_files, read.csv))
    final_file <- paste0("./new_data_yearly/", export_name, "lu_direct_results.csv")
    write.csv(all_results, final_file, row.names = FALSE)
    file.remove(unlist(csv_files)) # delete intermediate files
}


export_predicted_df <- function(iter) {
    i <- iterations$interval[iter]
    j <- iterations$data_par[iter]
    k <- iterations$biome[iter]
    l <- iterations$basic_par[iter]

    data <- dataframes[[i]][[k]]
    pars_names <- data_pars_names[[j]]
    biome_name <- biomes[[k]]
    basic_pars_name <- basic_pars_names[[l]]
    pars_iter <- find_combination_pars(iter, data)

    # Perform cross-validation and process results
    cv_output <- cross_valid(data, pars_iter, conditions)

    data[["pred"]] <- cv_output$pred
    data[["pred_nolag"]] <- cv_output$pred_nolag
    print(mean(data[["pred"]]))
    print(mean(data[["pred_nolag"]]))
    write.csv(data, paste0(c(import_name, biome_name, pars_names, "results.csv"), collapse = "_"), row.names = FALSE)
}


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


# Main

direct_optimization(iterations)
