
library(ggplot2)
library(foreach)
library(doParallel)
library(tidyverse)
library(mgcv)
library(randomForest)


# Source external R scripts for data import and function definitions
source("fit_1_import_data.r")
source("2_functions.r")

set.seed(1)
ncores <- 30
registerDoParallel(cores = ncores)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------------------- Switches ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fit_logistic <- FALSE

run_find_combination_pars <- TRUE

optim_switch <- TRUE
lm_switch <- FALSE
rf_switch <- FALSE
all_switch <- FALSE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Global Variables ------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# number of rows to be included in analysis
n_samples <- 10000


non_data_pars <- c("k0", "B0_exp", "B0", "theta")

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
lulc_categories <- list("aggregated", "non_aggregated")
combinations <- expand.grid(intervals, lulc_categories)
names <- apply(combinations, 1, paste, collapse = "_")

datafiles <- paste0("./data/", names, ".csv")
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparing dataframes whether or not there is a split between biomes

# 1 = Amazonia
# 2 = Caatinga
# 3 = Cerrado
# 4 = Mata Atlantica
# 5 = Pampa
# 6 = Pantanal
regions <- c("amaz", "atla", "both")

dataframes <- lapply(dataframes, function(df) {
    df[df$biome %in% c(1, 4), ]
})

# Split each dataframe by biome and store in a nested list
dataframes <- lapply(dataframes, function(df) {
    split(df, df$biome)
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
    "age|agbd|latitude|longitude|prec|si|mature_biomass|distance|biome|cwd",
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

if (fit_logistic) {
    basic_pars <- basic_pars[1:3]
}

basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))

# for linear model and random forest
# - remove climatic_pars (since they are only used in optim for yearly change)
# - add mature_biomass (nearest neighbor) to each parameter set
# - add new parameter sets
data_pars_lm <- c(
    lapply(
        Filter(function(x) !any(climatic_pars %in% x), data_pars), # remove climatic_pars
        function(x) c(x, "mature_biomass") # , "age")
    ),
    list(
        c("age"), c("mature_biomass"), c("age", "mature_biomass")
    )
)

filtered_data_pars_names <- data_pars_names[!sapply(data_pars, function(x) any(climatic_pars %in% x))]
data_pars_lm_names <- c(filtered_data_pars_names, "age", "mat_biomass", "age_mat_biomass")



pars_fit <- readRDS("./data/amaz_ideal_par_combination.rds")


data <- dataframes[[1]]
lu_pars <- pars_fit[[81]][!names(pars_fit[[81]]) %in% c("theta", "B0_exp", "k0")]
pars <- c(B0 = mean(dataframes[[1]][["agbd"]]), theta = 1, lu_pars) # age = 0, cwd = 0)

r_squared_values_opt <- c()
r_squared_values_lm <- c()

indices <- sample(c(1:5), nrow(data), replace = TRUE)

for (i in 1:5) {

    print(i)
    # Define the test and train sets
    test_data <- data[indices == i, ]
    train_data <- data[!indices == i, ]

    modellm <- run_lm(train_data, test_data, names(lu_pars))
    print(modellm$rsq)
    print(modellm$model_par)

    modelopt <- run_optim(train_data, test_data, pars, conditions)
    print(modelopt$rsq)
    print(modelopt$model_par)

}
 

