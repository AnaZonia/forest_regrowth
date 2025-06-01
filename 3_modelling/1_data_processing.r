# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                 Forest Regrowth Model Data Processing Functions
#
#                            Ana Avila - May 2025
#
#     This script defines the core functions used in the data processing and
#     preparation stages of the forest regrowth modeling process.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyverse)
library(fastDummies)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Data Import & Preparation Function ------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Function to import data, remove unnecessary columns, filter by biome,
# and optionally convert categorical variables to dummy variables.
#
# Arguments:
#   path             : The path to the CSV file.
#   biome            : Integer. Specifies which biome to retain (1 = Amazonia, 2 = Caatinga, etc.).
#   n_samples        : Integer. The maximum number of rows to sample from the filtered dataframe (default: 10,000).
#
# Returns:
#   df              : A dataframe filtered for the selected biome, cleaned, and optionally transformed with dummy variables.

import_data <- function(path, biome, n_samples = 10000, asymptote = "nearest_mature") {

    # if the folder exists, read all csv files in the folder and bind them together.
    # if the folder does not exist, read the single csv file.
    if (dir.exists(paste0("./0_data/", path))) {
        csv_files <- list.files(paste0("./0_data/", path), pattern = "\\.csv$", full.names = TRUE)
        df <- csv_files %>%
            map(read_csv) %>%
            bind_rows()
    } else {
        df <- read_csv(paste0("./0_data/", path, ".csv"))
        df <- df %>%
            select(-c(".geo", "system:index"))
    }

    # Convert categorical to factors
    df <- df %>%
        mutate(across(all_of(categorical), as.factor)) %>%
        filter(biome == biome)

    df <- na.omit(df)

    # remove columns with less than 100 non-zero values
    non_zero_counts <- colSums(df != 0, na.rm = TRUE)
    df <- df[, non_zero_counts > 100]

    # remove columns with less than 50 unique values
    df <- df %>%
        group_by(across(all_of(categorical))) %>%
        filter(n() >= 50) %>%
        ungroup() %>%
        mutate(across(all_of(categorical), droplevels))
    
    df <- dummy_cols(df,
        select_columns = categorical,
        remove_first_dummy = TRUE,
        remove_selected_columns = TRUE
    )

    asymptotes <- c("nearest_mature", "ecoreg_biomass", "quarter_biomass")

    if (asymptote == "full_amazon") {
        df$asymptote <- mean(df$nearest_mature, na.rm = TRUE)
        df <- df %>% select(-all_of(c(asymptotes, "quarter", "biome")))
    } else {
        # remove the columns in asymptotes that are not the designated asymptote
        # rename to asymptote the column named the same as the value of asymptote
        df <- df %>% rename(asymptote = !!sym(asymptote))
        df <- df %>% select(-all_of(c(
            asymptotes[asymptotes != asymptote],
            "quarter", "biome"
        )))
    }

    if (n_samples == "all") {
        coords <- df[, c("lat", "lon")]
        features <- df[, !names(df) %in% c("lat", "lon")]
        return(list(df = features, coords = coords))
    } else {

        if ("secondary_area" %in% names(df)) df <- df[, names(df) != "secondary_area"]

        df <- df[, !names(df) %in% c("lat", "lon")]
        df <- df[sample(nrow(df), min(n_samples, nrow(df)), replace = FALSE), ]
        return(df)
    }
}
