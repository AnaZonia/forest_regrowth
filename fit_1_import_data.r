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
    # Select numeric columns for normalization, excluding specified ones
    norm_cols <- c(names(data)[!grepl(paste0(c(unlist(categorical), "agbd", "nearest_mature_biomass"), collapse = "|"), names(data))])

    data <- data %>%
        mutate(across(
            all_of(norm_cols),
            # Normalization formula: (x - min(x)) / (max(x) - min(x))
            ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
        )) %>%
        # Remove columns that are entirely NA
        select(where(~ sum(is.na(.)) < nrow(data)))

    return(data)
}

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

import_data <- function(path, convert_to_dummy, process_climatic = TRUE) {
    
    columns_to_remove <- c(
        ".geo", "latitude", "longitude", "mature_forest_years",
        "last_LU",
        "num_fires_before_first_anthro", "num_fires_after_first_anthro", "num_fires_during_anthro"
    )

    df <- read_csv(path, show_col_types = FALSE) %>% # show_col_types = FALSE quiets a large message during import
        select(-any_of(columns_to_remove)) %>%
        select(-starts_with("system")) %>%
        mutate(across(all_of(categorical), as.factor))

    df <- df[df$biome %in% c(1, 4), ]
    
    if ("nearest_mature" %in% names(df)){
        names(df)[names(df) == "nearest_mature"] <- "nearest_mature_biomass"
    }

    list_of_dfs <- split(df, df$biome)
 
    list_of_dfs <- list(
        list_of_dfs[[1]], # First split result (e.g., biome 1)
        list_of_dfs[[2]], # Second split result (e.g., biome 4)
        df
    )

    list_of_dfs <- lapply(list_of_dfs, function(df) {

        if (process_climatic) {
            df <- process_climatic(df)
        }

        # df <- normalize(df)
        
        non_zero_counts <- colSums(df != 0, na.rm = TRUE)
        df <- df[, non_zero_counts > 50]

        df <- df %>%
            group_by(across(all_of(categorical))) %>%
            filter(n() >= 50) %>%
            ungroup() %>%
            mutate(across(all_of(categorical), droplevels))
        
        df <- df[sample(nrow(df), n_samples, replace = FALSE), ]
        df <- df[, !names(df) %in% "biome"]
        
        # Create dummy variables
        if (convert_to_dummy) {
            df <- dummy_cols(df,
                select_columns = categorical,
                remove_first_dummy = TRUE,
                remove_selected_columns = TRUE
            )
        }
        return(df)
    })


    print("Imported!")
    return(list_of_dfs)
}