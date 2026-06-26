# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#              Parameter Sets and Configurations
#
#                  Ana Avila - August 2025
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------------ Global Variables ----------------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

land_use <- c("lu", "fallow")
fires <- c("num_fires")
landscape <- c("dist", "sur_cover")
categorical <- c("topography", "last_lu", "ecoreg")
binary <- c("floodable_forests", "protec", "indig")
soil <- c("nitro", "phh2o", "ocd", "cec", "sand", "clay", "soc", "cfvo")

non_data_pars <- c("k0", "lag", "theta")

# Conditions for parameter constraints
conditions <- list('pars["k0"] < 0')

excluded_columns <- c("age", "biomass", "asymptote", "ecoreg", "area", "edge", "sd", "lat", "lon", "nearest_mature", "quarter", "biome")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Configurations of parameters ----------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

basic_pars_options <- list(
    lag = c("lag", "k0"),
    intercept = c("k0")
)

data_pars_options <- function(colnames) {
    return(list(
        age_only = c(),

        land_use = colnames[grepl(paste0(c(land_use, "last_lu"), collapse = "|"), colnames)],

        fires = colnames[grepl("num_fires", colnames)],

        environment = colnames[!grepl(paste0(c(excluded_columns, land_use, landscape, "num_fires"), collapse = "|"), colnames)],

        all = colnames[!grepl(paste0(c(excluded_columns), collapse = "|"), colnames)]
    ))
}

