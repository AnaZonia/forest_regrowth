
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
# categorical <- c("ecoreg", "soil", "indig", "protec")

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
intervals <- list("all")

datafiles <- paste0("./new_data/", name_import, "_", intervals, ".csv")
# datafiles <- paste0("./new_data/eu.csv")

dataframes <- lapply(datafiles, import_data, convert_to_dummy = TRUE, process_climatic = FALSE)

# dataframes[[1]][[1]] <- dataframes[[1]][[1]][, !grepl("ts|si|after|prec", names(dataframes[[1]][[1]]))]
# dataframes[[1]][[1]] <- dataframes[[1]][[1]][, !grepl("sur_cover", names(dataframes[[1]][[1]]))]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# --------------------------------- Define Parameters -----------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
land_use <- c("lulc", "LU", "fallow", "num_fires_before_regrowth", "distance", "sur_cover")

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

    biome_pars <- list(
        c(colnames_filtered[!grepl(paste0(categorical, collapse = "|"), colnames_filtered)]),# "cwd"),
        # c(colnames_filtered[!grepl(paste0(land_use, collapse = "|"), colnames_filtered)]),
        c(colnames_filtered) #, "cwd"),
    )

    data_pars[[i]] <- biome_pars
}


data_pars_names <- c(
    "no_categorical",
    # "no_land_use",
    "all"
)


# Define basic parameter sets for modeling
basic_pars <- list(
    c("age", "B0"),
    c("age", "k0", "B0"),
    c("k0", "B0"),
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

if (export_results) {
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


source("fit_2_functions.r")
install.packages("factoextra")
library(factoextra)


iter = 4
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


# Sample 10k unique rows for the first data frame
sample1_indices <- sample(nrow(data), size = 10000, replace = FALSE)
df1 <- data[sample1_indices, ]

# Remove the sampled rows from the original data, then sample another 10k rows for the second data frame
remaining_data <- data[-sample1_indices, ]
data2 <- sample(nrow(remaining_data), size = 10000, replace = FALSE)
data2 <- remaining_data[data2, ]
data <- df1
  
norm_data <- normalize_independently(data, data)$train_data
data_kmeans <- norm_data[, (names(norm_data) %in% names(pars_iter))]
kmeans_result <- kmeans(data_kmeans, 100, iter.max = 100)
# fviz_cluster(kmeans_result, data = data_kmeans)
fviz_nbclust(data_kmeans, kmeans, method = "silhouette")
# fviz_nbclust(data_kmeans, kmeans, method = "wss") +
#     geom_vline(xintercept = 3, linetype = 2)
# fviz_nbclust(data_kmeans, kmeans, nstart = 25, method = "gap_stat", nboot = 50)
norm_data$cluster <- kmeans_result$cluster

# Growth curve function
growth_curve <- function(pars, data) {
    B0 <- pars[["B0"]]
    k <- pars[["k"]]
    theta <- pars[["theta"]]

    return(B0 + (data[["nearest_mature_biomass"]] - B0) * (1 - exp(-k * data[["age"]]))^theta)
}


# Create a function to calculate residuals (difference between predicted and observed values)
residual_sum_of_squares <- function(pars, data, conditions) {
    predicted_agbd <- growth_curve(pars, data)
    residuals <- predicted_agbd - data[["agbd"]]

    result <- sum(residuals^2)

    # Check whether any of the parameters is breaking the conditions (e.g. negative values)
    if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}

# Optimization function
fit_growth_curve <- function(cluster_data) {
    # Initialize parameters
    initial_params <- c(
        "B0" = mean(data[["agbd"]]),
        "theta" = 1,
        "k" = 0.1
    ) # Initial guess for k

    conditions <- c(conditions, list('pars["k"] < 0'))
    conditions <- c(conditions, list('pars["B0"] < 0'))

    optim_result <- optim(
        par = initial_params,
        fn = residual_sum_of_squares,
        data = cluster_data,
        conditions = conditions
    )

    # Return the optimized parameters
    return(optim_result$par)
}

# Initialize a data frame to store results
k_values <- data.frame(cluster = numeric(), k = numeric())

for (c in 1:length(unique(norm_data$cluster))) {
    df_cluster <- subset(norm_data, cluster == c)
    # print(head(df_cluster))
    par <- fit_growth_curve(df_cluster)
    k_values <- rbind(k_values, data.frame(cluster = c, k = par["k"]))
    print(par["k"])
    pred <- growth_curve(par, df_cluster)
    print(calc_rsq(df_cluster, pred))
}

par <- fit_growth_curve(norm_data)
pred <- growth_curve(par, norm_data)
print(calc_rsq(norm_data, pred))


# Merge the k values into the original data by the cluster column
norm_data <- merge(norm_data, k_values, by = "cluster", all.x = TRUE)

# require(xgboost)
# library(Matrix)

sparse_data_kmeans <- Matrix(data.matrix(data_kmeans), sparse = TRUE)

bstSparse <- xgboost(data = sparse_data_kmeans, label = norm_data$k, max.depth = 2, eta = 1, nthread = 2, nrounds = 10)

new_k <- predict(bstSparse, sparse_data_kmeans)
norm_data$k <- new_k

# Growth curve function
growth_curve_2 <- function(pars, data) {
    B0 <- pars[["B0"]]
    k <- data[["k"]]
    theta <- pars[["theta"]]

    return(B0 + (data[["nearest_mature_biomass"]] - B0) * (1 - exp(-k * data[["age"]]))^theta)
}

# Create a function to calculate residuals (difference between predicted and observed values)
residual_sum_of_squares <- function(pars, data, conditions) {
    predicted_agbd <- growth_curve_2(pars, data)
    residuals <- predicted_agbd - data[["agbd"]]

    result <- sum(residuals^2)

    # Check whether any of the parameters is breaking the conditions (e.g. negative values)
    if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}

# Optimization function
fit_growth_curve <- function(cluster_data) {
    # Initialize parameters
    initial_params <- c(
        "B0" = mean(data[["agbd"]]),
        "theta" = 1)

    # conditions <- c(conditions, list('pars["k"] < 0'))
    conditions <- c(conditions, list('pars["B0"] < 0'))

    optim_result <- optim(
        par = initial_params,
        fn = residual_sum_of_squares,
        data = cluster_data,
        conditions = conditions
    )

    # Return the optimized parameters
    return(optim_result$par)
}

# Initialize a data frame to store results
k_values <- data.frame(cluster = numeric(), k = numeric())

norm_data2 <- normalize_independently(data2, data2)$train_data
data_kmeans <- norm_data2[, (names(norm_data2) %in% names(pars_iter))]
sparse_data_kmeans <- Matrix(data.matrix(data_kmeans), sparse = TRUE)
new_k <- predict(bstSparse, sparse_data_kmeans)
norm_data2$k <- new_k

par <- fit_growth_curve(norm_data2)
pred <- growth_curve_2(par, norm_data2)
print(calc_rsq(norm_data2, pred))






# Create a function to calculate residuals (difference between predicted and observed values)
residual_sum_of_squares <- function(pars, data, conditions) {
    predicted_agbd <- growth_curve(pars, data)
    residuals <- predicted_agbd - data[["agbd"]]

    result <- sum(residuals^2)

    # Check whether any of the parameters is breaking the conditions (e.g. negative values)
    if (any(sapply(conditions, function(cond) eval(parse(text = cond))))) {
        return(-Inf)
    } else if (is.na(result) || result == 0) {
        return(-Inf)
    } else {
        return(result)
    }
}
