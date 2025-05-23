# parameters.R - Define all parameter sets and configurations

# Global Variables
climatic_pars <- c("srad", "soil", "temp", "vpd", "aet", "def", "pr", "pdsi", "si")
# climatic_pars <- c("sdsr", "nsat", "musc")
land_use <- c("lu", "fallow", "num_fires")
landscape <- c("dist", "sur_cover")
categorical <- c("ecoreg", "topography") # , "last_lu")
binary <- c("floodable_forests", "protec", "indig")
non_data_pars <- c("k0", "B0", "lag", "theta")
interval <- c("5yr", "10yr", "15yr", "all")

# check permutation of
permut_check <- c("ecoreg", "num_fires", "dist", "sur_cover", "mean_srad", "mean_temp", "mean_aet")

# Conditions for parameter constraints
conditions <- list('pars["k0"] < 0')

excluded_columns <- c("age", "biomass", "nearest_mature")

# Configuration definition
basic_pars_options <- list(
    lag = c("lag", "k0"),
    intercept = c("B0", "k0")
)

data_pars_options <- function(colnames) {
    return(list(
        age_only = c(),

        # land_use_landscape = colnames[grepl(paste0(c(land_use, landscape), collapse = "|"), colnames)],

        environment = colnames[!grepl(paste0(c(excluded_columns, land_use, landscape, paste0(climatic_pars, "_")), collapse = "|"), colnames)],

        # environment_yearly = colnames[!grepl(paste0(c(excluded_columns, land_use, landscape, paste0("mean_", climatic_pars)), collapse = "|"), colnames)],

        # all_no_categorical = colnames[!grepl(paste0(c(excluded_columns, paste0(climatic_pars, "_"), categorical), collapse = "|"), colnames)], # all parameters, and no categorical variables

        all_mean_climate = colnames[!grepl(paste0(c(excluded_columns, paste0(climatic_pars, "_")), collapse = "|"), colnames)] # all parameters, and climatic variables as historical summaries,

        # all_yearly_climate = colnames[!grepl(paste(c(excluded_columns, paste0("mean_", climatic_pars)), collapse = "|"), colnames)]
        
        # all parameters, and each year's climatic variables included
        
        # all_yearly_climate = colnames[!grepl(paste0(c(excluded_columns, paste0("mean_", climatic_pars)), collapse = "|"), colnames)] # all parameters, and each year's climatic variables included
    ))
}
