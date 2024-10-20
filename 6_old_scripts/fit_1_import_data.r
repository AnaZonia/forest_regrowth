library(tidyverse)
library(fastDummies)

# - Imports the dataframe
# - Removes unnecessary columns that will not be used in analysis
# - Converts categorical data to dummy variables

#------------------ Global Variables ------------------#

climatic_pars <- c("prec", "si")

#------------------- Main Functions -------------------#

import_data <- function(path, convert_dummies = TRUE) {
    data <- read_csv(path, show_col_types = FALSE) #show_col_types = FALSE quiets a large message during import
    data <- subset(data, biome == 4)
    # data <- data %>%
    #     rename(nearest_mature_biomass = nearest_mature)

    data <- data[, !names(data) %in% c(
        "biome", "ecoreg", "topography", "cwd", "mean_aet"
    )]

    data <- data %>%
        mutate(across(
            where(is.numeric) &
                !matches("soil|biome|ecoreg|last_LU|protec|indig|agbd|nearest_mature_biomass"),
            ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
        )) %>%
        select(where(~ sum(is.na(.)) < nrow(data)))

    # Convert specified variables to factors
    categorical <- c("last_LU")
    data <- data %>%
        mutate(across(all_of(categorical), as.factor))
    data <- dummy_cols(data, select_columns = categorical, remove_selected_columns = TRUE)

    data <- sample_n(data, 10000)

    data
}


import_climatic_data <- function(data) {

    means <- sapply(climatic_pars, function(var) {
        rowMeans(data[, grep(var, names(data))],
            na.rm = TRUE
        )
    })

    colnames(means) <- paste0("mean_", climatic_pars)
    data <- cbind(data, means)

    df_climatic_hist <- tibble()
    for (age in 1:max(data$age)) {
        age_data <- data %>% filter(age == .env$age)
        years <- seq(2019, 2019 - age + 1, by = -1)
        # Identify all columns including the variables in climatic_pars
        clim_columns <- expand.grid(climatic_pars, years) %>%
            unite(col = "col", sep = "_") %>%
            pull(col)

        # subsect the dataframe to only include the climatic columns of the desired years
        all_clim_columns <- names(data)[str_detect(names(data), paste(climatic_pars, "_", collapse = "|"))]

        # turn all values in the columns of years not included to 0
        clim_columns_not_included <- setdiff(all_clim_columns, clim_columns)
        age_data[clim_columns_not_included] <- 0

        df_climatic_hist <- bind_rows(df_climatic_hist, age_data)
    }

    return(df_climatic_hist)

}
