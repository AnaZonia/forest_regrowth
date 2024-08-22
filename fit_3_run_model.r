
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

find_ideal_combination_pars <- TRUE

optim_switch <- FALSE
lm_switch <- FALSE
rf_switch <- FALSE
all_switch <- FALSE

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

# # Loop through the lists and dataframes to write CSV files
# for (i in seq_along(dataframes)) {
#     list_name <- intervals[i]
#     list_of_dfs <- dataframes[[i]]

#     for (j in seq_along(list_of_dfs)) {
#         df_name <- biomes[j]
#         df <- list_of_dfs[[j]]

#         # Create the filename
#         filename <- paste0(list_name, "_", df_name, ".csv")

#         # Write the CSV file
#         write.csv(df, filename, row.names = FALSE)

#         cat("Wrote file:", filename, "\n")
#     }
# }

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

if (fit_logistic) {
    basic_pars <- basic_pars[1:3]
}

basic_pars_names <- as.list(sapply(basic_pars, function(par_set) paste(par_set, collapse = "_")))

# for linear model and random forest
# - remove climatic_pars (since they are only used in optim for yearly change)
# - add nearest_mature (nearest neighbor) to each parameter set
# - add new parameter sets
data_pars_lm <- c(
    lapply(
        Filter(function(x) !any(climatic_pars %in% x), data_pars), # remove climatic_pars
        function(x) c(x, "nearest_mature") # , "age")
    ),
    list(
        c("age"), c("nearest_mature"), c("age", "nearest_mature")
    )
)

filtered_data_pars_names <- data_pars_names[!sapply(data_pars, function(x) any(climatic_pars %in% x))]
data_pars_lm_names <- c(filtered_data_pars_names, "age", "mat_biomass", "age_mat_biomass")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create grids of different combinations of model inputs

# Optim / growth curve
iterations_optim <- expand.grid(
    interval = seq_along(intervals),
    data_par = seq_along(data_pars),
    basic_par = seq_along(basic_pars),
    biome = seq_along(biomes)
)

basic_pars_with_age <- which(sapply(basic_pars, function(x) "age" %in% x))
data_pars_with_climatic <- which(sapply(data_pars, function(x) any(climatic_pars %in% x)))
# Remove rows where both conditions are met
iterations_optim <- iterations_optim %>% filter(!(basic_par %in% basic_pars_with_age & data_par %in% data_pars_with_climatic))

# Linear Model
iterations_lm <- expand.grid(
    interval = seq_along(intervals),
    data_par = seq_along(data_pars_lm)
)


# Random Forest
# Filter out single-predictor cases
rows_to_remove <- sapply(iterations_lm$data_par, function(i) {
    length(data_pars_lm[[i]]) == 1
})
iterations_rf <- iterations_lm[!rows_to_remove, ]

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





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------------------------- Test Area ---------------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

datafiles <- paste0("./data/amaz_", intervals, ".csv")
dataframes <- lapply(datafiles, import_climatic_data, normalize = TRUE)


pars_fit <- readRDS("./data/amaz_ideal_par_combination.rds")

data <- dataframes[[1]]
data <- data %>% rename(nearest_mature = mature_biomass)

for(i in seq_along(pars_fit)){
    lu_pars <- pars_fit[[i]][!names(pars_fit[[i]]) %in% non_data_pars]

    # pars <- c(theta = 1, lu_pars)
    # if ('B0' %in% names(pars_fit[[i]])){
    #     pars <- append(pars, B0 = mean(dataframes[[1]][["agbd"]]))
    # }
    # if ('B0' %in% names(pars_fit[[i]])){

    r_squared_values_opt <- c()
    r_squared_values_lm <- c()

    indices <- sample(c(1:5), nrow(data), replace = TRUE)

    for (i in 1:5) {
        print(i)
        # Define the test and train sets
        test_data <- data[indices == i, ]
        train_data <- data[!indices == i, ]

        modellm <- run_lm(train_data, names(lu_pars), test_data)
        print(modellm$rsq)
        print(modellm$model_par)
        r_squared_values_lm <- append(r_squared_values_lm, modellm$rsq)
        process_row(modellm)

        modelopt <- run_optim(train_data, pars_fit[[i]], test_data, conditions)
        print(modelopt$rsq)
        print(modelopt$model_par)
        r_squared_values_opt <- append(r_squared_values_opt, modelopt$rsq)
    }

}



pars <- pars_fit[[1]]
data = train_data
growth_curve(pars, data)
head(pars[["B0"]] + (data[["nearest_mature"]] - pars[["B0"]]) * (1 - exp(-k))^pars[["theta"]])
