# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#           Import and Normalization Functions
#
#                 Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(tidyverse)
library(fastDummies)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Data Import & Preparation Function ----------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' This function imports CSV files, cleans unnecessary columns, filters by biome, and handles categorical variables (via factor conversion and dummy encoding).
#'
#' @param path Character. Path to the folder containing CSV files.
#' @param biome_num Integer. Biome to retain (1 = Amazonia, 4 = Atlantic Forest).
#' @param n_samples Integer or "all". Maximum number of rows to sample
#'        (default = 10000). Use "all" to keep the entire dataset.
#' @param asymptote Character. Which asymptote column to use. Options:
#'        "nearest_mature", "ecoreg_biomass", "quarter_biomass", or "full_amazon".
#'
#' @return A data frame ready for modeling, or if \code{n_samples = "all"},
#'         a list with:
#'         \item{df}{Features (excluding coordinates).}
#'         \item{coords}{Latitude and longitude.}
#'
#' @details
#' - Removes columns with very few non-zero values or insufficient unique categories.
#' - Encodes categorical variables into dummy columns.
#' - Handles asymptote selection by renaming and filtering columns accordingly.
#' - Drops rows with missing values.
#'


import_data <- function(path, biome, n_samples = 10000, asymptote = "nearest_mature") {

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


    # if there is a column named distance_deep_forest rename it to dist
    if ("distance_deep_forest" %in% names(df)) {
        df <- df %>% rename(dist = distance_deep_forest)
    }

    # remove columns with less than 50 non-zero values
    df <- df %>% select(where(~ sum(. != 0) >= 50))

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Prepare Dataframes Function -----------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' Normalize numeric columns independently with min-max scaling
#'
#' Normalizes numerical columns of a dataset to the [0,1] range, while
#' excluding specified variables. If test data is provided, it is normalized
#' using the statistics from the training data.
#'
#' @param train_data Data frame. Training dataset.
#' @param test_data Optional data frame. Test dataset to normalize using the
#'        training data statistics.
#'
#' @return A list containing:
#' \item{train_data}{Normalized training dataset.}
#' \item{test_data}{(If provided) Normalized test dataset.}
#' \item{train_stats}{Data frame of min and max values used for scaling.}
#'

normalize_independently <- function(train_data, test_data = NULL) {
    # Identify numeric columns to normalize (excluding those in exclusion_list)
    exclusion_list <- c(categorical, binary, "biomass", "asymptote", "age", "area", "edge")

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
