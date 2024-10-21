
library(foreach)
library(doParallel)
library(tidyverse)
library(mgcv)
library(randomForest)


# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

set.seed(1)
ncores <- 50
registerDoParallel(cores = ncores)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Switches ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

find_ideal_combination_pars <- TRUE
export_results <- TRUE
old_data <- FALSE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

name_import <- "non_aggregated"
name_export <- name_import

# number of rows to be included in analysis
n_samples <- 10000
re_base <- dnorm(1000)

# List of climatic parameters that change yearly
climatic_pars <- c("prec", "si")
categorical <- c("ecoreg", "topography", "indig", "protec")
# categorical <- c("ecoreg", "soil", "last_LU", "indig", "protec")

# Define conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0')
non_data_pars <- c("k0", "B0_exp", "B0", "theta", "m_base", "sd_base")

# biomes <- c("amaz", "atla", "pant", "all")
# biomes <- c("amaz", "atla", "both")
biomes <- c("atla")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define land-use history intervals to import four dataframes
# intervals <- list("5yr", "10yr", "15yr", "all")
intervals <- list("15yr")

datafiles <- paste0("./new_data/", name_import, "_", intervals, ".csv")
# datafiles <- "./new_data/data_mapbiomas.csv"
source("fit_1_import_data.r")
dataframes <- lapply(datafiles, import_data, convert_to_dummy = TRUE, process_climatic = FALSE)

# dataframes_lm <- lapply(datafiles, import_data, convert_to_dummy = FALSE)

# dataframe lengths differ, all slightly below 10k. maybe check if that is making a difference in the results.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


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
        "age|agbd|prec|si|nearest_mature_biomass|cwd",
        colnames_intersect
    )]
    print(colnames_filtered)

    if (old_data) {
        biome_pars <- list(
            c("cwd", "mean_prec", "mean_si"),
            c(colnames_filtered[!grepl(paste0(categorical, collapse = "|"), colnames_filtered)]),
            c(colnames_filtered[!grepl(paste0(categorical, collapse = "|"), colnames_filtered)], "cwd", "mean_prec", "mean_si"),
            c(colnames_filtered, "cwd", "mean_prec", "mean_si")
        )
    } else {
        biome_pars <- list(
            c(colnames_filtered[!grepl(paste0(categorical, collapse = "|"), colnames_filtered)]),# "cwd"),
            c(colnames_filtered) #, "cwd")
        )
    }

    data_pars[[i]] <- biome_pars
}

if (old_data) {
    data_pars_names <- c(
        "mean_clim",
        "land_use",
        "land_use_clim",
        "land_use_clim_ecoreg_soil"
    )
} else {
    data_pars_names <- c(
        "no_categorical",
        "all"
    )
}

# Define basic parameter sets for modeling
basic_pars <- list(
    c("age", "B0"),
    c("age", "k0", "B0"),
    c("k0", "B0"),
    c("m_base", "sd_base", "k0")
)

basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))


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
# ----------------------------- Find ideal parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if (find_ideal_combination_pars) {
    # keep combination of parameters for initial fit
    initial_pars <- find_combination_pars(iterations_optim)
} else {
    initial_pars <- readRDS(paste0("./data/", name_export, "_ideal_par_combination.rds"))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source("fit_2_functions.r")


if (export_results) {
    # results <- foreach(
    #     iter = 1:nrow(iterations_optim),
    #     .combine = "bind_rows", .packages = c("dplyr", "randomForest")
    # ) %dopar% {
    for (iter in 1:length(initial_pars)){

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
        optim_cv_output <- cross_valid(data, run_optim, pars_iter, conditions)

        row <- process_row(optim_cv_output, "optim", intervals[[i]], pars_names, biome_name, basic_pars_name)
        print(row)
        row
    }

    write.csv(results, paste0("./data/", name_export, "_results.csv"), row.names = FALSE)
}

