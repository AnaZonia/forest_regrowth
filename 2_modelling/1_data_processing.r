# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#                 Forest Regrowth Model Data Processing Functions
#
#                            Ana Avila - May 2025
#
#     This script defines the core functions used in the data processing and
#     preparation stages of the forest regrowth modeling process.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyverse)
library(fastDummies)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# -------------------------- Data Import & Preparation Function ------------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

import_data <- function(path, biome_num, n_samples = 10000, asymptote = "nearest_mature") {

    csv_files <- list.files(paste0("./0_data/", path), pattern = "\\.csv$", full.names = TRUE)
    df <- csv_files %>%
        map(read_csv) %>%
        bind_rows()

    # Convert categorical to factors
    df <- df %>%
        mutate(across(any_of(categorical), as.factor)) %>%
        filter(biome == biome_num)

    # remove columns with less than 100 non-zero values
    non_zero_counts <- colSums(df != 0, na.rm = TRUE)
    df <- df[, non_zero_counts > 100]

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

    # remove any rows with NA values (important due to gaps in CMIP6 data)
    df <- df %>% filter(rowSums(is.na(.)) == 0)

    # if there is a column named distance_deep_forest rename it to dist
    if ("distance_deep_forest" %in% names(df)) {
        df <- df %>% rename(dist = distance_deep_forest)
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Prepare Dataframes Function -----------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Function used in import_data to normalize numeric columns in dataframes.
# Arguments:
#   train_data  : The training dataframe to be used for analysis
#   test_data   : (Optional) The testing dataframe to be normalized using the same
#                 parameters as train_data.
# Returns:
#   A list containing:
#     - train_data : A dataframe with normalized numerical values
#     - test_data  : (If provided) The test dataframe, normalized using train_data statistics

normalize_independently <- function(train_data, test_data = NULL) {
    # Identify numeric columns to normalize (excluding those in exclusion_list)
    exclusion_list <- c(categorical, binary, "biomass", "asymptote", "age")

    norm_cols <- names(train_data)[!grepl(paste0(exclusion_list, collapse = "|"), names(train_data))]

    # Compute summary statistics
    train_stats <- data.frame(
        variable = norm_cols,
        min = sapply(norm_cols, function(var) min(train_data[[var]], na.rm = TRUE)),
        max = sapply(norm_cols, function(var) max(train_data[[var]], na.rm = TRUE))
    )

    # Apply Min-Max scaling using the precomputed min and max
    for (i in seq_along(norm_cols)) {
        var <- norm_cols[i]
        train_data[[var]] <- (train_data[[var]] - train_stats$min[i]) /
            (train_stats$max[i] - train_stats$min[i])

        if (!is.null(test_data)) {
            test_data[[var]] <- (test_data[[var]] - train_stats$min[i]) /
                (train_stats$max[i] - train_stats$min[i])
        }
    }

    if (is.null(test_data)) {
        return(list(train_data = train_data, train_stats = train_stats))
    } else {
        # keep in test_data only rows with values greater than zero (those with values in the range of the training data)
        return(list(train_data = train_data, test_data = test_data, train_stats = train_stats))
    }
}
