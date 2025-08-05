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
data_10k <- import_data("grid_10k_amazon_secondary", biome = 1, n_samples = 20000, asymptote = "nearest_mature")

norm_data <- normalize_independently(data_10k)
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

data_1k <- import_data("grid_1k_amazon_secondary", biome = 1, n_samples = 20000, asymptote = "nearest_mature")
norm_data <- normalize_independently(data_1k)
norm_data <- norm_data$train_data

basic_pars <- basic_pars_options[["lag"]]
data_pars <- c()

init_pars <- find_combination_pars(basic_pars, data_pars, norm_data)
model <- run_optim(norm_data, init_pars, conditions)
pred <- growth_curve(model$par, norm_data, lag = model$par["lag"])
r2 <- calc_r2(norm_data, pred)
print(r2)




tst_df_clean <- function(df){
    # Convert categorical to factors
    df <- df %>%
        mutate(across(any_of(categorical), as.factor)) %>%
        filter(biome == 1)

    # # remove columns with less than 100 non-zero values
    # non_zero_counts <- colSums(df != 0, na.rm = TRUE)
    # df <- df[, non_zero_counts > 100]

    # # remove columns with less than 50 unique values
    # df <- df %>%
    #     group_by(across(any_of(categorical))) %>%
    #     # filter(n() >= 50) %>%
    #     ungroup() %>%
    #     mutate(across(any_of(categorical), droplevels))

    # df <- dummy_cols(df,
    #     select_columns = categorical,
    #     remove_first_dummy = TRUE,
    #     remove_selected_columns = TRUE
    # )

    # asymptotes <- c("nearest_mature", "ecoreg_biomass", "quarter_biomass")

    # if (asymptote == "full_amazon") {
    #     df$asymptote <- mean(df$nearest_mature, na.rm = TRUE)
    #     df <- df %>% select(-any_of(c(asymptotes, "quarter", "biome")))
    # } else {
    #     # remove the columns in asymptotes that are not the designated asymptote
    #     # rename to asymptote the column named the same as the value of asymptote
    #     df <- df %>% rename(asymptote = !!sym(asymptote))
    #     df <- df %>% select(-any_of(c(
    #         asymptotes[asymptotes != asymptote],
    #         "quarter", "biome"
    #     )))
    # }

    # # remove any rows with NA values (important due to gaps in CMIP6 data)
    # df <- df %>% filter(rowSums(is.na(.)) == 0)
    return(df)
}

csv_files <- list.files(paste0("./0_data/grid_1k_amazon_secondary"), pattern = "\\.csv$", full.names = TRUE)
data_1k <- csv_files %>%
    map(read_csv) %>%
    bind_rows()

csv_files <- list.files(paste0("./0_data/grid_10k_amazon_secondary"), pattern = "\\.csv$", full.names = TRUE)
data_10k <- csv_files %>%
    map(read_csv) %>%
    bind_rows()

data_1k_tst <- tst_df_clean(data_1k)
data_10k_tst <- tst_df_clean(data_10k)

mean(data_10k_tst$biomass, na.rm = TRUE)
mean(data_1k_tst$biomass, na.rm = TRUE)

sd(data_10k$biomass, na.rm = TRUE)
sd(data_1k$biomass, na.rm = TRUE)

nrow(data_10k_tst)
nrow(data_1k_tst)

# randomly select from data_1k the same number of rows as in data_10k
set.seed(1)
data_1k_sampled <- data_1k_tst %>%
    sample_n(nrow(data_10k_tst))

nrow(data_1k_sampled)
mean(data_1k_sampled$biomass, na.rm = TRUE)
