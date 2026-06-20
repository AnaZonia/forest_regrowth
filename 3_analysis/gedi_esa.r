# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#      Run comparisons for answering to reviews
#
#                     Ana Avila - May 2026
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(foreach)
library(doParallel)
library(tidyverse)
library(xtable)

source("2_modelling/1_parameters.r")
source("2_modelling/1_data_processing.r")
source("2_modelling/2_modelling.r")
source("2_modelling/2_cross_validate.r")
source("2_modelling/2_forward_selection.r")

# Set up parallel processing
set.seed(1)
ncore <- 4
registerDoParallel(cores = ncore)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------- GEDI / ESA CCI comparison ---------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



import_data <- function(path, biome, n_samples = 10000) {
    csv_files <- list.files(paste0("./0_data/", path), pattern = "\\.csv$", full.names = TRUE)

    df <- csv_files %>%
        map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
        bind_rows()
    
    # remove columns ESA_biomass, nearest_mature
    # df <- df %>% select(-ESA_biomass, -nearest_mature, -biome)
    # df <- df %>% rename(asymptote = "nearest_mature_GEDI")

    df <- df %>% select(-biomass, -nearest_mature_GEDI, -biome)
    df <- df %>% rename(asymptote = "nearest_mature", biomass = "ESA_biomass")

    # remove columns with all NA values
    df <- df[, colSums(is.na(df)) < nrow(df)]
    # remove any rows with NA values
    df <- df %>% filter(rowSums(is.na(.)) == 0)

    # Convert categorical to factors
    df <- df %>%
        mutate(across(any_of(categorical), as.factor)) %>%
        filter(biome == biome)

    # remove columns with less than 50 unique values
    df <- df %>%
        group_by(across(any_of(categorical))) %>%
        ungroup() %>%
        mutate(across(any_of(categorical), droplevels))

    df <- dummy_cols(df,
        select_columns = categorical,
        remove_first_dummy = TRUE,
        remove_selected_columns = TRUE
    )


    # remove columns with less than 50 non-zero values
    df <- df %>% select(where(~ sum(. != 0) >= 50))

    df <- subset(df, df$biomass < 600)

    if (n_samples == "all") {
        coords <- df[, c("lat", "lon")]
        features <- df[, !names(df) %in% c("lat", "lon")]
        return(list(df = features, coords = coords))
    } else {
        if ("area" %in% names(df)) df <- df[, names(df) != "area"]
        df <- df[, !names(df) %in% c("lat", "lon")]
        df <- df[sample(nrow(df), min(n_samples, nrow(df)), replace = FALSE), ]
        return(df)
    }
}



data <- import_data("grid_10k_secondary_GEDI_ESA_edges_removed", biome = 1, n_samples = 30000)

basic_pars <- basic_pars_options[["lag"]]
data_pars <- data_pars_options(colnames(data))[["all"]]

cv_results <- cross_validate(data, basic_pars, data_pars, conditions, 5)

cv_results

cv_results_gedi <- cv_results



# remove from data the rows where biomass is greater than 600
# data <- subset(data, data$biomass < 600)

names(df)

# tst <- read.csv("./0_data/gedi_esa/age_gedi_esa.csv")

summary(lm(biomass ~ age, data = df))

summary(lm(ESA_biomass ~ age, data = df))

# result is the same with the data extracted for the same locations.



