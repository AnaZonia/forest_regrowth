
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#      Compare R2 of different asymptotes and land use aggregations
#
#                     Ana Avila - August 2025
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
# --------- One lag per ecoregion of the Amazon ----------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


asymptote <- "nearest_mature"
path = "uncertainty_propagation"

csv_files <- list.files(paste0("./0_data/", path), pattern = "\\.csv$", full.names = TRUE)

df <- csv_files %>%
    map(~ suppressMessages(read_csv(.x, show_col_types = FALSE, progress = FALSE))) %>%
    bind_rows()

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
    select_columns = "topography",
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
)

asymptotes <- c("nearest_mature", "ecoreg_biomass", "quarter_biomass")

if (asymptote == "full_amazon") {
    df$asymptote <- mean(df$nearest_mature, na.rm = TRUE)
    df <- df %>% select(-any_of(c(asymptotes, "quarter", "biome")))
} else {
    # remove the columns in asymptotes that are not the designated asymptote
    # rename to asymptote the column named the same as the value of asymptote
    df <- df %>% rename(asymptote = !!sym(asymptote))
    df <- df %>% select(-any_of(c(
        asymptotes[asymptotes != asymptote],
        "quarter", "biome"
    )))
}

# remove columns with less than 50 non-zero values
df <- df %>% select(where(~ sum(. != 0) >= 50))


df <- subset(df, ave(seq_along(ecoreg), ecoreg, FUN = length) >= 500)

# names(df)

df <- df[, -which(names(df) %in% c("lat", "lon", "area", "sd"))]

results <- c()

for (ecoregion in unique(df$ecoreg)) {
    print(ecoregion)

    df_ecoreg <- subset(df, df$ecoreg == ecoregion)

    basic_pars <- basic_pars_options[["lag"]]
    data_pars <- data_pars_options(colnames(df_ecoreg))[["all"]]

    cv_results <- cross_validate(df_ecoreg, basic_pars, data_pars, conditions, folds = 2)

    print(cv_results[[3]])

    results <- c(results, mean(cv_results[[3]]))
}




