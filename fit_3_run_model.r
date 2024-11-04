
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
ncores <- 5
registerDoParallel(cores = ncores)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

name_import <- "non_aggregated"
name_export <- name_import

# number of rows to be included in analysis
n_samples <- 10000
re_base <- 0 # rnorm(10)

# List of climatic parameters that change yearly
climatic_pars <- c("srad", "soil", "temp", "vpd")
categorical <- c("ecoreg", "topography")

# Define conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0', 'pars["k0"] < 0')
non_data_pars <- c("k0", "B0", "theta", "m_base", "sd_base")

# biomes <- c("amaz", "atla", "both")
biomes <- c("atla")
# biomes <- c("amaz")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define land-use history intervals to import four dataframes
# intervals <- list("5yr", "10yr", "15yr", "all")
intervals <- list("all")

datafiles <- paste0("./new_data_yearly/", name_import, "_", intervals, ".csv")
dataframes <- lapply(datafiles, import_data, convert_to_dummy = TRUE)
dataframes_lm <- lapply(datafiles, import_data, convert_to_dummy = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


land_use <- c("lulc", "LU", "fallow", "num_fires")
landscape <- c("distance", "sur_cover", "nearest_mature_biomass")

# Identify common columns across all dataframes of the same biome (across intervals)
# (avoids errors for parametesrs like last_LU, that not all dataframes may contain depending on how many rows were sampled)
# and filter out non-predictors
unique_colnames <- c()
data_pars <- c()
for (i in seq_along(biomes)){
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
        c(colnames_filtered)
    )

    # biome_pars <- list(
    #     c(colnames_filtered[!grepl(paste0(c(land_use, landscape), collapse = "|"), colnames_filtered)]),
    #     c(colnames_filtered[!grepl(paste0(land_use, collapse = "|"), colnames_filtered)]),
    #     c(colnames_filtered[!grepl(paste0(landscape, collapse = "|"), colnames_filtered)]),
    #     c(colnames_filtered),
    #     c(colnames_filtered_no_mean_climate, climatic_pars)
    # )

    data_pars[[i]] <- biome_pars
}

data_pars_names <- c(
    "all"
)

# data_pars_names <- c(
#     "no_land_use_no_landscape",
#     "no_land_use",
#     "no_landscape",
#     "all",
#     "all_yearly_clim"
# )

# Define basic parameter sets for modeling
basic_pars <- list(
    # c("age", "k0", "B0"), #age is fit as its own independent predictor
    # c("k0", "B0"), # k is multiplied by the age column
    c("m_base", "sd_base", "k0")
)

basic_pars_names <- c(
    # "intercept_age",
    # "intercept",
    "lag"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

iterations <- expand.grid(
    basic_par = seq_along(basic_pars),
    interval = seq_along(intervals),
    data_par = seq_along(data_pars[[1]]),
    biome = seq_along(biomes)
)

basic_pars_fit_age <- which(sapply(basic_pars, function(x) "age" %in% x))
data_pars_with_climatic <- which(sapply(data_pars[[1]], function(x) any(climatic_pars %in% x)))
# Remove rows where both conditions are met
iterations <- iterations %>% filter(!(basic_par %in% basic_pars_fit_age & data_par %in% data_pars_with_climatic))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# source("fit_2_functions.r")

pars_iter_list <- list()
for (iter in 1:nrow(iterations)) {
    i <- iterations$interval[iter]
    k <- iterations$biome[iter]
    data <- dataframes[[i]][[k]]
    pars_iter <- find_combination_pars(iter, data)
    pars_iter_list[[iter]] <- pars_iter
    saveRDS(pars_iter_list, "temp_pars_iter_list.rds")
}

# Define batch size
batch_size <- 20
num_batches <- ceiling(nrow(iterations) / batch_size)

# Initialize an empty list to store file names
csv_files <- vector("list", num_batches)

# Loop over each batch
for (batch in 1:num_batches) {
    start <- (batch - 1) * batch_size + 1
    end <- min(batch * batch_size, nrow(iterations))
    print(start)
    print(end)

    # Process the current batch
    results <- foreach(iter = start:end, .combine = rbind, .packages = c("tidyverse")
    ) %dopar% {
    # for (iter in start:end){

        i <- iterations$interval[iter]
        j <- iterations$data_par[iter]
        k <- iterations$biome[iter]
        l <- iterations$basic_par[iter]

        pars_iter <- pars_iter_list[[iter]]
        data <- dataframes[[i]][[k]]
        pars_names <- data_pars_names[[j]]
        biome_name <- biomes[[k]]
        basic_pars_name <- basic_pars_names[[l]]

        # Perform cross-validation and process results
        cv_output <- cross_valid(data, pars_iter, conditions)

        row <- process_row(cv_output, basic_pars_name, intervals[[i]], pars_names, biome_name)

        if (!any(climatic_pars %in% names(pars_iter))) {
            data_lm <- dataframes_lm[[i]][[k]]
            # Perform cross-validation and process results
            cv_output <- cross_valid(data_lm, pars_iter)
            row_lm <- process_row(cv_output, "lm", intervals[[i]], pars_names, biome_name)
            row <- rbind(row, row_lm)
        }
        print(row)
        return(row)
    }

    # Define the name for the intermediate CSV file
    batch_file <- paste0("./new_data_yearly/", name_export, "_batch_", batch, "_results.csv")
    csv_files[[batch]] <- batch_file

    # Write the batch results to a CSV file
    write.csv(results, batch_file, row.names = FALSE)
}

# Combine all intermediate CSV files and export final CSV
all_results <- do.call(rbind, lapply(csv_files, read.csv))
final_file <- paste0("./new_data_yearly/", name_export, "_final_results.csv")
write.csv(all_results, final_file, row.names = FALSE)



# lag, atla, all, all
n_clusters <- 5

for (iter in 1:nrow(iterations)){
    # Extract iteration-specific parameters
    iter = 1
    i <- iterations$interval[iter]
    j <- iterations$data_par[iter]
    k <- iterations$biome[iter]
    l <- iterations$basic_par[iter]

    data <- dataframes[[i]][[k]]
    data_pars_name <- data_pars_names[[j]]
    biome_name <- biomes[[k]]
    basic_pars_name <- basic_pars_names[[l]]

    data$pred <- NA
    # names_anthro <- names(data)[grepl(paste0(c(land_use, landscape), collapse = "|"), names(data))]
    names_enviro <- names(data)[!grepl(paste0(c(land_use, landscape), collapse = "|"), names(data))]
    data_norm <- data[, names(data) %in% names_enviro]
    data_norm <- normalize_independently(NULL, data_norm)$train_data
    pca <- PCA(data_norm, ncp = n_clusters, graph = FALSE)
  
    data_norm_pca <- as.data.frame(pca$ind$coord)
    kmeans_result <- kmeans(data_norm_pca, centers = n_clusters, nstart = 20)
    data$cluster <- kmeans_result$cluster

    r2_cluster <- numeric(n_clusters)
    # Loop through each unique cluster and perform cross-validation

    for (cluster_id in 1:n_clusters) {
        print(cluster_id)
        data_cluster <- subset(data, cluster == cluster_id)    
        pars_iter <- find_combination_pars(iter, data_cluster)
        # Perform cross-validation and process results
        optim_cv_output <- cross_valid(data_cluster, pars_iter, conditions)
        print(optim_cv_output$r2_mean)
        print(optim_cv_output$r2_sd)
        r2_cluster[cluster_id] <- optim_cv_output$r2_final
        data[data$cluster == cluster_id, c("pred_cv", "pred_final")] <- optim_cv_output$pred
        print(mean(optim_cv_output$pred))
    }

    print(mean(r2_cluster))
    
    head(data[, c("pred_cv", "pred_final")])
    mean(data[["pred_cv"]])
    mean(data[["pred_final"]])
    # row <- cbind(row, data.frame(r2_mean_clusters = mean(r2_cluster)))
    # print(row)
    # return(row)
}


