# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                 Forest Regrowth Model Data Processing Functions
#
#                            Ana Avila - March 2025
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
#   convert_to_dummy : Logical (TRUE/FALSE). If TRUE, categorical variables are converted to dummy variables.
#   columns_to_remove : Vector of column names to exclude from the dataframe.
#   biome            : Integer. Specifies which biome to retain (1 = Amazonia, 2 = Caatinga, etc.).
#   n_samples        : Integer. The maximum number of rows to sample from the filtered dataframe (default: 10,000).
#
# Returns:
#   df              : A dataframe filtered for the selected biome, cleaned, and optionally transformed with dummy variables.

import_data <- function(path, convert_to_dummy, columns_to_remove, biome, n_samples = 10000) {
    
    df <- read_csv(datafiles[[1]], show_col_types = FALSE) %>% # show_col_types = FALSE quiets a large message during import
        mutate(across(all_of(categorical), as.factor)) %>%
        filter(biome == 4) %>%
        select(-all_of(c(columns_to_remove, "biome")))

    # remove columns with less than 100 non-zero values
    non_zero_counts <- colSums(df != 0, na.rm = TRUE)
    df <- df[, non_zero_counts > 100]

    # remove columns with less than 50 unique values
    df <- df %>%
        group_by(across(all_of(categorical))) %>%
        filter(n() >= 50) %>%
        ungroup() %>%
        mutate(across(all_of(categorical), droplevels))
    
    # Create dummy variables
    if (convert_to_dummy) {
        df <- dummy_cols(df,
            select_columns = categorical,
            remove_first_dummy = TRUE,
            remove_selected_columns = TRUE
        )
    }

    df <- df[sample(nrow(df), min(n_samples, nrow(df)), replace = FALSE), ]


    print("Imported!")
    return(df)
}