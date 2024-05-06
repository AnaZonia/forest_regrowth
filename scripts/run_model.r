####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
####################################################################
# Ana Avila - May 2024
# ~~~~~~~~~~~~~~~~~~~~
# Intakes:
# Outputs:
####################################################################

library(ggplot2)
library(terra)
library(tidyverse) # for stringr operations
library(mlr) # for createDummyFeatures
setwd("C:/Users/anaca/Desktop/forest_regrowth/")
source("./scripts/model_functions.r")

"
  - Imports the dataframe
  - Removes unnecessary columns that will not be used in analysis
  - Converts categorical data to dummy variables
"

years <- seq(1985, 2019, 1)

data <- read.csv('./data/final_w_biomass.csv')
# Drop unnecessary columns
data <- data[ , -which(names(data) %in% c("system.index",".geo", 'biome'))]
# create dummy variables for the categorical data with more than 2 types
categorical = c('ecoreg', 'soil')
data[categorical] <- lapply(data[categorical], as.factor)
data <- createDummyFeatures(data, cols = categorical)

# define the climatic parameters - the ones that change yearly
climatic <- c('prec', 'si')
# define the non-climatic parameters - the ones that are fixed throughout regrowth and
# that are used for fitting the model (excludes age and agbd)
non_climatic <- names(data)[!grepl("prec|si|agbd", names(data))]

# define the list that will hold the parameters and their values
pars <- c(B0 = 40, A = 80, theta = 1.5, sd = 1) # intercept, asymptote, shape term, standard deviation
pars_init <- c(climatic, non_climatic)
# define initial values for all parameters
pars_init_values <- rep(0.1, length(pars_init))
new_elements <- setNames(pars_init_values, pars_init)
pars <- c(pars, new_elements)

# write the growth curve with yearly climatic data and permanent non-climatic data
growth_curve <- function(pars, data) {
  k = data[[1]]*0
  for (year in years){
    for (clim_var in climatic){
      k <- k + pars[clim_var]*data[[paste0(clim_var, '_', year)]]
    }
    for (unique_var in non_climatic) {
      k <- k + pars[unique_var] * data[[unique_var]]
    }
  }
  
  return(pars['B0'] + pars['A'] * (1 - exp(-k))^pars['theta'])
}

o = optim(pars, fn = nls, data = data, G = growth_curve)
o
