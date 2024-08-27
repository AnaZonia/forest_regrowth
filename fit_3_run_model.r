
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

name <- "non_aggregated"

# number of rows to be included in analysis
n_samples <- 10000

# Define climatic parameters that change yearly
climatic_pars <- c("prec", "si")

# Define parameters that do not correspond to data, used for functional form
if (fit_logistic) {
    name <- paste0(name, "_logistic")
    conditions <- list()
    non_data_pars <- c("k0", "B0")
} else {
    # Define conditions for parameter constraints
    conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0')
    non_data_pars <- c("k0", "B0_exp", "B0", "theta")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ---------------------------------- Import Data ----------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define land-use history intervals to import four dataframes
intervals <- list("5yr", "10yr", "15yr", "all")
# lulc_categories <- list("aggregated", "non_aggregated")
# combinations <- expand.grid(intervals, lulc_categories)
# names <- apply(combinations, 1, paste, collapse = "_")

datafiles <- paste0("./data/", name, "_", intervals, ".csv")
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparing dataframes whether or not there is a split between biomes

# 1 = Amazonia
# 2 = Caatinga
# 3 = Cerrado
# 4 = Mata Atlantica
# 5 = Pampa
# 6 = Pantanal
biomes <- c("amaz", "atla", "both")

dataframes <- lapply(dataframes, function(df) {
    df[df$biome %in% c(1, 4), ]
})

# Split each dataframe by biome and store in a nested list
dataframes <- lapply(dataframes, function(df) {
    split_result <- split(df, df$biome)
    list(
        split_result[[1]], # First split result (biome 1)
        split_result[[2]], # Second split result (biome 4)
        df # Original dataframe
    )
})

dataframes <- lapply(dataframes, function(list_of_dfs) {
    lapply(list_of_dfs, function(df) {
        df_sampled <- df[sample(nrow(df), n_samples, replace = FALSE), ]
        return(df_sampled)
    })
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Identify common columns across all dataframes
# (avoids errors for parameters like last_LU, that not all dataframes may contain depending on how many rows were sampled)
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
    colnames_filtered[!grepl("ecoreg|soil", colnames_filtered)], # land use only
    colnames_filtered,
    c(colnames_filtered, "cwd", "mean_prec", "mean_si"),
    c(colnames_filtered, "cwd", climatic_pars)
)

data_pars_names <- c(
    "cwd", "meanclim", "yearlyclim", "LU", "LU_soil_ecor",
    "LU_soil_ecor_meanclim", "LU_soil_ecor_yearlyclim"
)

# Define basic parameter sets for modeling
basic_pars <- list(
    c("age", "B0"),
    c("age", "k0", "B0"),
    c("k0", "B0"), # age is implicit
    c("k0", "B0_exp") # age is implicit
)
basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))

if (fit_logistic) {
    basic_pars <- basic_pars[1:3]
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

# Optim / growth curve
iterations_optim <- expand.grid(
    interval = seq_along(intervals),
    data_par = seq_along(data_pars),
    biome = seq_along(biomes),
    basic_par = seq_along(basic_pars)
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
    initial_pars <- readRDS(paste0("./data/", name, "_ideal_par_combination.rds"))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Run Model ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

results <- foreach(
    iter = 1:nrow(iterations_optim),
    .combine = "bind_rows", .packages = c("dplyr", "randomForest")
    ) %dopar% {
    # Extract iteration-specific parameters
    i <- iterations_optim$interval[iter]
    j <- iterations_optim$data_par[iter]
    k <- iterations_optim$biome[iter]
    l <- iterations_optim$basic_par[iter]
    
    data <- dataframes[[i]][[k]]
    pars_names <- data_pars_names[[j]]
    biome_name <- biomes[[k]]
    basic_pars_names <- basic_pars_names[[l]]
    pars_iter <- initial_pars[[iter]] # Parameters obtained from "find combination pars"

    # Perform cross-validation and process results
    cross_valid <- cross_valid(data, pars_iter, conditions)
    row <- process_row(cross_valid[[1]], "optim", intervals[[i]], pars_names, biome_name,
        basic_pars_names = basic_pars_names
    )

    if (!any(climatic_pars %in% names(pars_iter))) {

        row_lm <- process_row(cross_valid[[2]], "lm", intervals[[i]], pars_names, biome_name)
        row <- rbind(row, row_lm)

    #     if (length(lu_pars) > 1) {
    #         # Perform cross-validation and process results
    #         cross_valid_rf <- cross_valid(data, run_rf, unique(c(lu_pars, "nearest_mature")))
    #         run_rf <- process_row(cross_valid_rf, "rf", intervals[[i]], pars_names, biome_name)
    #         row <- rbind(row, run_rf)
    #     }
    }

    print(row)
    row
}

# write.csv(results, "lm_optim_results.csv", row.names = FALSE)


source("fit_2_functions.r")
data <- dataframes[[1]][[1]]
pars_iter <- initial_pars[[10]]

cross_valid1 <- cross_valid(data, pars_iter, conditions)
cross_valid1[[1]]
