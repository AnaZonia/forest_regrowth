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
        ".geo", "latitude", "longitude", "mature_forest_years", "fallow", "last_LU",
        "num_fires_before_first_anthro", "num_fires_after_first_anthro", "num_fires_during_anthro"
    )

    data <- read_csv(path, show_col_types = FALSE) %>% # show_col_types = FALSE quiets a large message during import
        select(-any_of(columns_to_remove)) %>%
        select(-starts_with("system")) %>%
        # select(-starts_with("LU")) %>%
        mutate(across(all_of(categorical), as.factor))

    data <- data[data$biome %in% c(1, 4), ]

    list_of_dfs <- split(data, data$biome)

    list_of_dfs <- list(
        list_of_dfs[[1]], # First split result (e.g., biome 1)
        list_of_dfs[[2]], # Second split result (e.g., biome 4)
        data # Original dataframe without split
    )

    list_of_dfs <- lapply(list_of_dfs, function(df) {
        df_sampled <- df[sample(nrow(df), n_samples, replace = FALSE), ]
        df_sampled <- df_sampled %>%
            group_by(across(all_of(categorical))) %>% # Group by the categorical columns
            filter(n() >= 200) %>% # Keep groups with 50 or more occurrences
            ungroup() %>% # Ungroup to return to a normal dataframe
            mutate(across(all_of(categorical), droplevels))
        df_sampled <- process_climatic(df_sampled)
        df_sampled <- normalize(df_sampled)

        # Create dummy variables
        if (convert_to_dummy) {
            df_sampled <- dummy_cols(df_sampled, select_columns = categorical, remove_selected_columns = TRUE)
        }
        return(df_sampled)
    })

    print("Imported!")
    return(list_of_dfs)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Process Climatic Variables ---------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Function to import and process climatic data, normalize it, and optionally convert categorical variables to dummy variables.
#
# Arguments:
#   data   : A list of dataframes to process.
# Returns:
#   df_climatic_hist : A processed dataframe with historical climatic data.

process_climatic <- function(data) {

    # Calculate the mean of specified climatic variables across years
    means <- sapply(climatic_pars, function(var) {
        rowMeans(data[, grep(var, names(data))],
            na.rm = TRUE
        )
    })
    # Append mean climatic variables as new columns to the data
    colnames(means) <- paste0("mean_", climatic_pars)
    data <- cbind(data, means)

    # Process climatic data for each age group
    df_climatic_hist <- tibble()


    for (yrs in 1:max(data$age)) {
        age_data <- data %>% filter(age == yrs)

        # Generate a sequence of years for the current age group
        # Starting from 2019 and going back 'yrs' number of years
        year_seq <- seq(2019, 2019 - yrs + 1, by = -1)

        # Create column names for climatic parameters for each year
        clim_columns <- expand.grid(climatic_pars, year_seq) %>%
            unite(col = "col", sep = "_") %>%
            pull(col)

        # Subset the dataframe to only include the climatic columns of the desired years
        all_clim_columns <- names(data)[str_detect(names(data), paste(climatic_pars, "_", collapse = "|"))]

        # Identify climatic columns not relevant for the current age group
        clim_columns_not_included <- setdiff(all_clim_columns, clim_columns)
        # Set values in non-relevant columns to 0
        age_data[clim_columns_not_included] <- 0

        df_climatic_hist <- bind_rows(df_climatic_hist, age_data)
    }

    return(df_climatic_hist)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Prepare Dataframes Function --------------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Function used in import_data to normalize numeric columns in dataframes.
# Arguments:
#   data             : The dataframe to be used for analysis
# Returns:
#   data             : A dataframe with normalized numerical values


normalize <- function(data) {
    data <- data %>%
        mutate(across(
            where(is.numeric) &
                !matches("soil|biome|ecoreg|last_LU|protec|indig|agbd|nearest_mature|fallow"),
            # Normalization formula: (x - min(x)) / (max(x) - min(x))
            ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
        )) %>%
        # Remove columns that are entirely NA
        select(where(~ sum(is.na(.)) < nrow(data)))

    return(data)
}

