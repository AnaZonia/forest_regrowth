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
categorical <- c("ecoreg", "topography", "last_lu")
binary <- c("floodable_forests", "protec", "indig")
soil <- c("nitro", "phh2o", "ocd", "cec", "sand", "clay", "soc", "cfvo")

non_data_pars <- c("k0", "B0", "lag")

# Conditions for parameter constraints
conditions <- list('pars["k0"] < 0')

# "mean_def", "mean_temp", "mean_pr", "phh2o" excluded due to multicollinearity
excluded_columns <- c("age", "biomass", "asymptote", "mean_def", "mean_temp", "mean_pr", "phh2o")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ------------ Configurations of parameters ----------------#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

basic_pars_options <- list(
    lag = c("lag", "k0"),
    intercept = c("B0", "k0")
)

data_pars_options <- function(colnames) {
    return(list(
        age_only = c(),

        land_use = colnames[grepl(paste0(c(land_use, "last_lu"), collapse = "|"), colnames)],

        fires = colnames[grepl("num_fires", colnames)],

        environment = colnames[!grepl(paste0(c(excluded_columns, land_use, landscape, "num_fires"), collapse = "|"), colnames)],

        all_mean_climate = colnames[!grepl(paste0(c(excluded_columns), collapse = "|"), colnames)]
    ))
}

