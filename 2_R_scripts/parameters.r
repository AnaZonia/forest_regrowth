# parameters.R - Define all parameter sets and configurations

# Global Variables
climatic_pars <- c("srad", "soil", "temp", "vpd", "aet", "def", "pdsi", "pr")
land_use <- c("lu", "fallow", "num_fires")
landscape <- c("dist", "sur_cover")
categorical <- c("ecoreg", "topography") # , "last_lu")
non_data_pars <- c("k0", "B0", "theta", "lag")
columns_to_remove <- c("ecoreg_biomass", "quarter_biomass", "quarter", ".geo", "system.index")

# Conditions for parameter constraints
conditions <- list('pars["theta"] > 10', 'pars["theta"] < 0', 'pars["k0"] < 0')

# Configuration definition
get_config <- function() {
    list(
        basic_pars_options = list(
            lag = c("lag", "k0", "theta"),
            intercept = c("k0", "B0", "theta"),
            intercept_age_multiplicative = c("k0", "B0", "theta", "age")
        ),
        data_pars_options = c("land_use_landscape_only", "all"),
        biomes = 1:3, # If you have multiple biomes
        ncores = 4,
        n_samples = 15000
    )
}

# Generate data parameter sets based on names and data columns
generate_data_pars <- function(data_pars_name, colnames) {
    exclusion_pattern <- paste(c(
        "age", "biomass", "nearest_biomass",
        paste0(climatic_pars, "_")
    ), collapse = "|")

    if (data_pars_name == "land_use_only") {
        data_pars <- colnames[grepl(paste0(c(land_use), collapse = "|"), colnames)]
    } else if ("landscape_only") { 
        data_pars <- landscape
    } else if (data_pars_name == "all") {
        data_pars <- colnames[!grepl(exclusion_pattern, colnames)]
    }

    return(data_pars)
}

# Get basic parameters based on selection
get_basic_pars <- function(basic_pars_name) {
    config <- get_config()
    return(config$basic_pars_options[[basic_pars_name]])
}