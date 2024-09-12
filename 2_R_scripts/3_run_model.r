
library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)
library(mgcv)
library(randomForest)


# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("fit_2_functions.r")

set.seed(1)
ncores <- 20
registerDoParallel(cores = ncores)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Switches ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fit_logistic <- FALSE
find_ideal_combination_pars <- FALSE
export_results <- FALSE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

name_import <- "non_aggregated"
name_export <- name_import

# number of rows to be included in analysis
n_samples <- 10000

# Define climatic parameters that change yearly
climatic_pars <- c("prec", "si")

# Define parameters that do not correspond to data, used for functional form
if (fit_logistic) {
    name_export <- paste0(name_import, "_logistic")
    conditions <- list()
    non_data_pars <- c("k0", "B0")
} else {
    # Define conditions for parameter constraints
    conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0')
    non_data_pars <- c("k0", "B0_exp", "B0", "theta")
}

biomes <- c("amaz", "atla", "both")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define land-use history intervals to import four dataframes
intervals <- list("5yr", "10yr", "15yr", "all")
# lulc_categories <- list("aggregated", "non_aggregated")
# combinations <- expand.grid(intervals, lulc_categories)
# names <- apply(combinations, 1, paste, collapse = "_")

datafiles <- paste0("./data/", name_import, "_", intervals, ".csv")
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE, convert_to_dummy = TRUE)
dataframes <- prepare_dataframes(dataframes, c(1, 4))
dataframes_lm <- lapply(datafiles, import_climatic_data, normalize = TRUE, convert_to_dummy = FALSE)
dataframes_lm <- prepare_dataframes(dataframes_lm, c(1, 4))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Identify common columns across all dataframes
# (avoids errors for parametesrs like last_LU, that not all 
# dataframes may contain depending on how many rows were sampled)
# and filter out non-predictors

colnames_lists <- lapply(dataframes, function(df_list) {
    lapply(df_list, function(df) colnames(df))
})

colnames_lists <- lapply(colnames_lists, unlist)

# Step 2: Find the intersection of all column names
colnames_intersect <- Reduce(intersect, colnames_lists)

colnames_filtered <- colnames_intersect[!grepl(
    "age|agbd|latitude|longitude|prec|si|nearest_mature|distance|biome|cwd",
    colnames_intersect
)]

# Define different sets of data parameters for modeling
data_pars <- list(
    c("cwd"),
    c("cwd", "mean_prec", "mean_si"),
    c("cwd", climatic_pars),
    c(colnames_filtered[!grepl("ecoreg|soil|last_LU", colnames_filtered)], "cwd", "mean_prec", "mean_si"), # land use only
    c(colnames_filtered, "cwd", "mean_prec", "mean_si"),
    c(colnames_filtered[!grepl("ecoreg|soil|last_LU", colnames_filtered)], "cwd", climatic_pars), # land use only
    c(colnames_filtered, "cwd", climatic_pars)
)

data_pars_names <- c(
    "cwd", "mean_clim", "yearly_clim", "all_continuous_mean_clim", "all_mean_clim",
    "all_continuous_yearly_clim", "all_yearly_clim"
)

if (!fit_logistic) {
    # Define basic parameter sets for modeling
    basic_pars <- list(
        c("age", "B0"),
        c("age", "k0", "B0"),
        c("k0", "B0"), # age is implicit
        c("k0", "B0_exp") # age is implicit
    )
    basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

iterations_optim <- expand.grid(
    interval = seq_along(intervals),
    data_par = seq_along(data_pars),
    biome = seq_along(biomes)
)

if (!fit_logistic) {
    basic_par_grid <- expand.grid(
        basic_par = seq_along(basic_pars)
    )
    iterations_optim <- merge(iterations_optim, basic_par_grid, by = NULL)
    
    basic_pars_with_age <- which(sapply(basic_pars, function(x) "age" %in% x))
    data_pars_with_climatic <- which(sapply(data_pars, function(x) any(climatic_pars %in% x)))
    # Remove rows where both conditions are met
    iterations_optim <- iterations_optim %>% filter(!(basic_par %in% basic_pars_with_age & data_par %in% data_pars_with_climatic))
}


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

if (export_results){
    results <- foreach(
        iter = 1:nrow(iterations_optim),
        .combine = "bind_rows", .packages = c("dplyr", "randomForest")
    ) %dopar% {
        # for(iter in 1:nrow(iterations_optim)){

        # Extract iteration-specific parameters
        i <- iterations_optim$interval[iter]
        j <- iterations_optim$data_par[iter]
        k <- iterations_optim$biome[iter]

        data <- dataframes[[i]][[k]]
        pars_names <- data_pars_names[[j]]
        biome_name <- biomes[[k]]

        if (!fit_logistic) {
            l <- iterations_optim$basic_par[iter]
            basic_pars_names <- basic_pars_names[[l]]
        }

        pars_iter <- initial_pars[[iter]] # Parameters obtained from "find combination pars"

        # Perform cross-validation and process results
        optim_cv_output <- cross_valid(data, run_optim, pars_iter, conditions)
        row <- process_row(optim_cv_output, "optim", intervals[[i]], pars_names, biome_name)

        if (!any(climatic_pars %in% names(pars_iter))) {
            data <- dataframes_lm[[i]][[k]]

            lu_pars <- names(pars_iter[!names(pars_iter) %in% non_data_pars])

            lm_cv_output <- cross_valid(data, run_lm, unique(c(lu_pars, "nearest_mature")))

            row_lm <- process_row(lm_cv_output, "lm", intervals[[i]], pars_names, biome_name)
            row <- rbind(row, row_lm)

        }

        print(row)
        row
    }

    write.csv(results, "optim_results.csv", row.names = FALSE)
}

