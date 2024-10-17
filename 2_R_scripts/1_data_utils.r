library(tidyverse)
library(fastDummies)

# - Imports the dataframe
# - Removes unnecessary columns that will not be used in analysis
# - Converts categorical data to dummy variables

#------------------ Global Variables ------------------#
# List of climatic parameters that change yearly
climatic_pars <- c("prec", "si")

#------------------- Main Functions -------------------#


# Function to import data and optionally convert categorical columns to dummy variables.
# Arguments:
#   path             : The path to the CSV file.
#   convert_to_dummy : Logical switch to decide if categorical variables should be converted to dummy variables.
#                      If TRUE, categorical variables will be split into dummy variables.
#                      If FALSE, categorical variables will remain as factors.
# Returns:
#   data             : A dataframe ready for analysis with or without dummy variables.

import_data <- function(path, convert_to_dummy = TRUE) {
    categorical <- c("ecoreg", "soil", "topography", "last_LU")

    # Import data
    data <- read_csv(path, show_col_types = FALSE)

    # Create dummy variables only for categorical columns that are present
    if (convert_to_dummy) {
        # Identify which categorical columns are present in the data
        present_categorical <- intersect(categorical, names(data))

        if (length(present_categorical) > 0) {
            # Convert present categorical columns to factors
            data <- data %>%
                mutate(across(all_of(present_categorical), as.factor))

            # Create dummy variables for present categorical columns
            data <- dummy_cols(data,
                select_columns = present_categorical,
                remove_selected_columns = TRUE
            )
        }
    }

    print("Imported!")
    return(data)
}



# Arguments:
#   path       : The path to the CSV file.
#   normalize  : Logical switch to decide if the data should be normalized.
#                If TRUE, numeric columns will be normalized.
#   convert_to_dummy : Logical switch passed to import_data function to handle categorical variables.
# Returns:
#   df_climatic_hist : A processed dataframe with historical climatic data.

import_climatic_data <- function(path, normalize, convert_to_dummy) {

    data <- import_data(path, convert_to_dummy)

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



# ------------------- Prepare Dataframes Function -------------------#
# This function prepares a list of dataframes by:
# 1. Filtering for specific biomes.
# 2. Splitting the dataframes by biome.
# 3. Sampling a specified number of rows from each split dataframe.
#

# 1 = Amazonia
# 2 = Caatinga
# 3 = Cerrado
# 4 = Mata Atlantica
# 5 = Pampa
# 6 = Pantanal

# Arguments:
#   df_list   : A list of dataframes to process.
#   filter_biomes : A vector of biome codes to filter on (e.g., c(1, 4) for Amazonia and Mata Atlantica).
#
# Returns:
#   A nested list of dataframes, split and sampled according to the specified biomes.

prepare_dataframes <- function(df_list, filter_biomes) {
    # Filter dataframes to include only rows matching the specified biomes
    filtered_dfs <- lapply(df_list, function(df) {
        df[df$biome %in% filter_biomes, ]
    })

    # Split each filtered dataframe by biome and store in a nested list
    split_dfs <- lapply(filtered_dfs, function(df) {
        split_result <- split(df, df$biome)
        list(
            split_result[[1]], # First split result (e.g., biome 1)
            split_result[[2]], # Second split result (e.g., biome 4)
            df # Original dataframe without split
        )
    })

    # Sample n_samples rows from each split dataframe
    sampled_dfs <- lapply(split_dfs, function(list_of_dfs) {
        lapply(list_of_dfs, function(df) {
            df_sampled <- df[sample(nrow(df), n_samples, replace = FALSE), ]
            # df_sampled[, colSums(df_sampled != 0) > 0]
            return(df_sampled)
        })
    })

    return(sampled_dfs)
}
