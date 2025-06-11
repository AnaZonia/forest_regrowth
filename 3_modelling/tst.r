# main.R - Main script
# Runs experiments and saves results for different parameter combinations

# Load required libraries
library(foreach)
library(doParallel)
library(tidyverse)

# Source other scripts
source("3_modelling/1_parameters.r")
source("3_modelling/1_data_processing.r")
source("3_modelling/2_modelling.r")
source("3_modelling/2_cross_validate.r")
source("3_modelling/2_feature_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# quick R2 check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 20000, asymptote = "nearest_mature")
norm_data <- normalize_independently(data)
norm_data <- norm_data$train_data

basic_pars <- basic_pars_options[["lag"]]
data_pars <- c()

init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
r2 <- calc_r2(norm_data, pred)
print(r2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# quick R2 check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("grid_1k_amazon_secondary", biome = 1, n_samples = 20000, asymptote = "nearest_mature")
norm_data <- normalize_independently(data)
norm_data <- norm_data$train_data

basic_pars <- basic_pars_options[["lag"]]
data_pars <- c()

init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
r2 <- calc_r2(norm_data, pred)
print(r2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# quick R2 check
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data <- import_data("grid_1k_amazon_secondary", biome = 1, n_samples = 20000, asymptote = "nearest_mature")
norm_data <- normalize_independently(data)
norm_data <- norm_data$train_data

basic_pars <- basic_pars_options[["lag"]]
data_pars <- c()

init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
r2 <- calc_r2(norm_data, pred)
print(r2)
# 0.3425
# 0.3559

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Mean and Yearly Climate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


data <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 10000, asymptote = "quarter_biomass")
norm_data <- normalize_independently(data)
norm_data <- norm_data$train_data

climatic_pars <- c("srad", "temp", "def", "vpd", "pr", "pdsi", "aet")
data_pars_options <- list(
    mean_climate = colnames(data)[grepl(paste(c(paste0("mean_", climatic_pars)), collapse = "|"), colnames(data))]
    # yearly_climate = climatic_pars
)


for (data_pars_name in names(data_pars_options)) {
    print(data_pars_name)
    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options[[data_pars_name]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
    r2 <- calc_r2(norm_data, pred)
    print(r2)
}

for (data_pars_name in names(data_pars_options)) {
    print(data_pars_name)
    basic_pars <- basic_pars_options[["intercept"]]
    data_pars <- data_pars_options[[data_pars_name]]
    init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
    model <- run_optim(norm_data, init_pars, conditions)
    pred <- growth_curve(model$par, norm_data)
    r2 <- calc_r2(norm_data, pred)
    print(r2)
}



data1k <- import_data("grid_1k_amazon_secondary", biome = 1, n_samples = 20000, asymptote = "nearest_mature")

# data1k <- list.files(paste0("./0_data/grid_1k_amazon_secondary"), pattern = "\\.csv$", full.names = TRUE) %>%
#     map(read_csv) %>%
#     bind_rows()

data10k <- import_data("grid_10k_amazon_secondary", biome_num = 4, n_samples = 10000, asymptote = "nearest_mature")

table(data10k$biome)


df <- list.files(paste0("./0_data/grid_10k_amazon_secondary"), pattern = "\\.csv$", full.names = TRUE) %>%
    map(read_csv) %>%
    bind_rows()

biome_num <- 4

# Convert categorical to factors
df2 <- df %>%
    mutate(across(any_of(categorical), as.factor)) %>%
    filter(biome == biome_num)

table(df$biome)

# print the column names in df that contain at least one NA value
na_cols <- sapply(df, function(x) any(is.na(x)))
na_cols <- names(df)[na_cols]
# Print columns with NAs
print(na_cols)

head(df$ecoreg_biomass)

# Find common columns
common_cols <- intersect(names(data1k), names(data10k))

# Subset to common columns and add group labels
data1k_sub <- data1k[, common_cols]
data10k_sub <- data10k[, common_cols]
data1k_sub$group <- "1k"
data10k_sub$group <- "10k"

# Combine
combined <- rbind(data1k_sub, data10k_sub)
