# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                 Forest Regrowth Model Data Processing Functions
#
#                            Ana Avila - August 2024
#
#     This script defines the core functions used in the data processing and
#     preparation stages of the forest regrowth modeling process.
#
#     Functions included:
#     - process_climatic
#     - normalize
#     - import_data
#
#     These functions handle data import, climatic variable processing,
#     normalization, and optional conversion of categorical variables to
#     dummy variables.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyverse)
library(fastDummies)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Prepare Dataframes Function --------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Function to import data, remove unnecessary columns, and optionally convert categorical columns to dummy variables.
# This function prepares a list of dataframes by:
# 1. Filtering for specific biomes.
# 2. Splitting the dataframes by biome.
# 3. Sampling a specified number of rows from each split dataframe.
#
# Arguments:
#   path             : The path to the CSV file.
#   convert_to_dummy : Logical switch to decide if categorical variables should be converted to dummy variables.
#                      If TRUE, categorical variables will be split into dummy variables.
#                      If FALSE, categorical variables will remain as factors.
# 1 = Amazonia
# 2 = Caatinga
# 3 = Cerrado
# 4 = Mata Atlantica
# 5 = Pampa
# 6 = Pantanal
#
# Returns:
#   list_of_dfs             : A list with three dataframes, one per ecoregion, ready for analysis with or without dummy variables.

import_data <- function(path, convert_to_dummy) {
    
    columns_to_remove <- c(
        ".geo", "latitude", "longitude", "pr_", "si_", "aet_", "last_LU", "GEDI"
    )

    df <- read_csv(datafiles[[1]], show_col_types = FALSE) %>% # show_col_types = FALSE quiets a large message during import
        select(-starts_with(columns_to_remove)) %>%
        mutate(across(all_of(categorical), as.factor))

    list_of_dfs <- split(df, df$biome)
 
    list_of_dfs <- list(
        list_of_dfs[[1]], # First split result (e.g., biome 1)
        list_of_dfs[[2]], # Second split result (e.g., biome 4)
        df
    )

    list_of_dfs <- lapply(list_of_dfs, function(df) {
        df <- df[, !names(df) %in% "biome"]
        non_zero_counts <- colSums(df != 0, na.rm = TRUE)
        df <- df[, non_zero_counts > 100]

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
        # df <- df[sample(nrow(df), min(n_samples, nrow(df)), replace = FALSE), ]

        return(df)
    })

    print("Imported!")
    return(list_of_dfs)
}