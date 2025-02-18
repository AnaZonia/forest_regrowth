
library(foreach)
library(doParallel)
library(tidyverse)
library(FactoMineR) # For PCA
library(factoextra) # Optional, for PCA visualization
library(cluster) # For k-means clustering

# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

set.seed(1)
ncores <- 20
registerDoParallel(cores = ncores)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Switches ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

find_ideal_combination_pars <- TRUE
export_results <- TRUE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

name_import <- "non_aggregated"
name_export <- name_import

# number of rows to be included in analysis
n_samples <- 10000
re_base <- 0 # rnorm(10)

# List of climatic parameters that change yearly
climatic_pars <- c("prec", "si")
categorical <- c("ecoreg", "topography", "indig", "protec", "last_LU")

# Define conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0')
non_data_pars <- c("k0", "B0_exp", "B0", "theta", "m_base", "sd_base")

# biomes <- c("amaz", "atla", "pant", "all")
biomes <- c("amaz", "atla")
# biomes <- c("atla")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("fit_1_import_data.r")

# Define land-use history intervals to import four dataframes
# intervals <- list("5yr", "10yr", "15yr", "all")
intervals <- list("all")

# datafiles <- paste0("./new_data/", name_import, "_", intervals, ".csv")
# dataframes <- lapply(datafiles, import_data, convert_to_dummy = TRUE, process_climatic = FALSE)

dataframe <- read.csv("~/Documents/data/mapbiomas_heinrich_field.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
land_use <- c("lulc", "LU", "fallow", "num_fires_before_regrowth")
landscape <- c("distance", "sur_cover", "nearest_mature_biomass")

# Identify common columns across all dataframes of the same biome (across intervals)
# (avoids errors for parametesrs like last_LU, that not all dataframes may contain depending on how many rows were sampled)
# and filter out non-predictors
data_pars <- c()
for (i in seq_along(biomes)){
    colnames_lists <- lapply(dataframes, function(df_list) {
        colnames(df_list[[i]])
    })

    # Step 2: Find the intersection of all column names
    colnames_intersect <- Reduce(intersect, colnames_lists)

    colnames_filtered <- colnames_intersect[!grepl(
        "age|agbd|prec|si|nearest_mature_biomass",
        colnames_intersect
    )]
    print(colnames_filtered)

    biome_pars <- list(
        c(colnames_filtered[!grepl(paste0(land_use, collapse = "|"), colnames_filtered)]),
        c(colnames_filtered[!grepl(paste0(landscape, collapse = "|"), colnames_filtered)]),
        c(colnames_filtered)
    )

    data_pars[[i]] <- biome_pars
}

data_pars_names <- c(
    "no_land_use",
    "no_landscape",
    "all"
)

# Define basic parameter sets for modeling
basic_pars <- list(
    c("age", "B0"),
    # c("age", "k0", "B0"),
    # c("k0", "B0"),
    c("m_base", "sd_base", "k0")
)

basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))
data_pars

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

iterations_optim <- expand.grid(
    basic_par = seq_along(basic_pars),
    interval = seq_along(intervals),
    data_par = seq_along(data_pars[[1]]),
    biome = seq_along(biomes)
)

basic_pars_with_age <- which(sapply(basic_pars, function(x) "age" %in% x))
data_pars_with_climatic <- which(sapply(data_pars, function(x) any(climatic_pars %in% x)))
# Remove rows where both conditions are met
iterations_optim <- iterations_optim %>% filter(!(basic_par %in% basic_pars_with_age & data_par %in% data_pars_with_climatic))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


if (not_cluster) {
    initial_pars <- find_combination_pars(iterations_optim)

    results <- foreach(
        iter = 1:nrow(iterations_optim),
        .combine = "bind_rows", .packages = c("dplyr", "randomForest")
    ) %dopar% {
        # for (iter in 1:length(initial_pars)){

        # Extract iteration-specific parameters
        i <- iterations_optim$interval[iter]
        j <- iterations_optim$data_par[iter]
        k <- iterations_optim$biome[iter]
        l <- iterations_optim$basic_par[iter]

        data <- dataframes[[i]][[k]]
        pars_names <- data_pars_names[[j]]
        biome_name <- biomes[[k]]
        basic_pars_name <- basic_pars_names[[l]]
        pars_iter <- initial_pars[[iter]] # Parameters obtained from "find combination pars"

        # Perform cross-validation and process results
        optim_cv_output <- cross_valid(data, pars_iter, conditions)


        row <- process_row(optim_cv_output, intervals[[i]], pars_names, biome_name, basic_pars_name)
        print(row)
        row
    }

    write.csv(results, paste0("./data/", name_export, "_results.csv"), row.names = FALSE)
}




n_clusters <- 10
if (export_results) {

    results <- data.frame(
        biome_name = character(),
        data = character(),
        data_pars = character(),
        basic_pars = character(),
        rsq_all = numeric(),
        rsq_mean = numeric(),
        rsq_sd = numeric(),
        stringsAsFactors = FALSE
    )

    for (iter in 1:nrow(iterations_optim)){

        # Extract iteration-specific parameters
        i <- iterations_optim$interval[iter]
        j <- iterations_optim$data_par[iter]
        k <- iterations_optim$biome[iter]
        l <- iterations_optim$basic_par[iter]

        data <- dataframes[[i]][[k]]
        data_pars_name <- data_pars_names[[j]]
        biome_name <- biomes[[k]]
        basic_pars_name <- basic_pars_names[[l]]

        data$pred <- NA
        data_norm <- data[, names(data) %in% data_pars[[1]][[j]]]
        data_norm <- normalize_independently(data_norm)$train_data
        pca <- PCA(data_norm, ncp = n_clusters, graph = FALSE)

        data_norm_pca <- as.data.frame(pca$ind$coord)
        kmeans_result <- kmeans(data_norm_pca, centers = n_clusters, nstart = 20)
        data$cluster <- kmeans_result$cluster

        cluster_r_squared <- numeric(n_clusters)
        # Loop through each unique cluster and perform cross-validation
        for (cluster_id in 1:n_clusters) {
            data_cluster <- subset(data, cluster == cluster_id)

            pars_iter <- find_combination_pars(iter, data_cluster)

            # Perform cross-validation and process results
            optim_cv_output <- cross_valid(data_cluster, pars_iter, conditions)
            print(optim_cv_output$rsq_mean)
            print(optim_cv_output$rsq_sd)
            cluster_r_squared[cluster_id] <- optim_cv_output$rsq_final
            data$pred[data$cluster == cluster_id] <- optim_cv_output$pred
        }

        row <- data.frame(
            biome_name = biome_name,
            data = intervals[[i]],
            data_pars = data_pars_name,
            basic_pars = basic_pars_name,
            rsq_mean = mean(cluster_r_squared),
            rsq_sd = sd(cluster_r_squared),
            stringsAsFactors = FALSE
        )

        # Append the row to results
        results <- bind_rows(results, row)

        print(row)
    }

    # write.csv(results, paste0("./data/", name_export, "_results.csv"), row.names = FALSE)
}

