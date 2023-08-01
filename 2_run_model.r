####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the district of Paragominas, in the Western Amazon.
# Ana Avila - Jan 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################

library(ggplot2)
library(terra) # handling spatial data
library(tidyverse)
library(sf)
setwd("/home/aavila/forest_regrowth")
source("/home/aavila/forest_regrowth/scripts/0_forest_regrowth_functions.r")

######################################################################
#################      unify data into central df   ##################
######################################################################
prec <- rast('./model_ready_rasters/prec_0000000000-0000095232.tif')
temp <- rast('./model_ready_rasters/temp_0000000000-0000095232.tif')
fire <- rast('./model_ready_rasters/0000000000-0000095232_fire_history_santoro.tif')
last_LU <- rast('last_LU_0000000000-0000095232_santoro.tif')
regrowth_paper <- rast('./secondary_forest_age_v2_2018/secondary_forest_age_v2_2018-0000000000-0000065536.tif')
regrowth <- crop(regrowth_paper, fire)
regrowth[regrowth == 0] <- NA

santoro_raster <- rast('./santoro/N00W050_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv4.0.tif')
santoro_raster <- crop(santoro_raster, regrowth)
santoro_raster <- resample(santoro_raster, regrowth,'near')
santoro_raster <- mask(santoro_raster, regrowth)
all_data <- c(santoro_raster, regrowth)

temp_cropped <- crop(temp, all_data)
prec_cropped <- crop(prec, all_data)

total_temp <- app(temp_cropped, sum) #what if I try to incorporate the yearly effect right now?
total_prec <- app(prec_cropped, sum)
names(fire_cropped) <- c('total_fires', 'ts_fire')
names(total_temp) <- c('total_temp')
names(total_prec) <- c('total_prec')

all_data_santoro <- c(all_data, total_prec, total_temp, fire, last_LU)
file_name <- paste0(tempfile(), "_.tif")
lu_tile <- makeTiles(all_data, c(2000,2000), file_name, na.rm = TRUE, overwrite=TRUE)
lu_tile <- lapply(lu_tile, rast)

rast_to_df <- function(raster){
  coords <- crds(raster, df=FALSE, na.rm=TRUE)
  values_stack <- terra::extract(raster, coords, cells=FALSE, method="simple")
  central_df <- values_stack[complete.cases(values_stack), ]
  return(central_df)
}

all_data_csv <- lapply(lu_tile, rast_to_df)
all_data_csv <- bind_rows(all_data_csv)
colnames(all_data_csv) <- c('agbd', 'age', 'prec', 'temp', 'total_fires', 'ts_fire', 'last_LU')

saveRDS(all_data_csv, 'santoro_ESA_alldata.rds')

agbdagb <- all_data_csv
all_data_csv <- readRDS('santoro_ESA_alldata.rds')

###############################################################

colnames(all_data_csv) <- c('agbd', 'age')
new_all_data <- aggregate(agbd~age, data=all_data_csv, median)
plot(new_all_data$age, new_all_data$agbd)

######################################################################
#################        Sampling     ##################
######################################################################
total_temp <- app(temp, sum)
total_prec <- app(prec, sum)

# getting percent of mature forest cover within x neighboring patches:
# would like to get values within 300m radius (range of dispersal - cite)
mature_mask <- rast('./mapbiomas/mature_masks/0000000000-0000095232_mature_mask.tif')
mature_mask <- terra::crop(mature_mask, regrowth)
mature_mask[mature_mask > 0] <- 1
mature_mask <- subst(mature_mask, NA, 0)

range <- 21
mature_sum <- focal(mature_mask, range, sum, na.rm = TRUE) # run focal on these areas
mature_sum <- mask(mature_sum, regrowth) # select only the sum of areas surrounding secondary patches.

# mean biomass of surrounding mature forests around 3km (~100 pixels of 30m)
biomass_range <- 101 # window size
mature_biomass <- mask(santoro_raster, mature_mask)  # get only biomass of mature forests
mature_total_biomass <- focal(mature_biomass, biomass_range, sum, na.rm = TRUE)


##############


# min-max normalize central_df
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

all_data_csv <- all_data_csv[all_data_csv$last_LU %in% c(15, 41, 48),]

all_data_csv$last_LU <- factor(all_data_csv$last_LU)
dummy_LU <- as.data.frame(model.matrix(~ all_data_csv$last_LU - 1))
names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')

all_data_csv <- all_data_csv[,-7]
normalized_df <- as.data.frame(lapply(all_data_csv, minMax))
data <- cbind(normalized_df, dummy_LU)

#####################################



fit <- lm(central_df$agbd ~ central_df$age)
summary(fit)
plot(central_df$agbd ~ central_df$age)
abline(fit)

sds <- aggregate(agbd ~ age, central_df, sd)
means <- aggregate(agbd ~ age, central_df, median)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'median', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = median)) + 
  geom_errorbar(aes(ymin = median - sd,
                    ymax = median + sd)) +
  geom_point() + theme(text = element_text(size = 20))  



######################################################################
#################        passing into the model     ##################
######################################################################
# create another list that sets to zero if things are not present

G <- function(pars) {
  # Extract parameters of the model
  E = pars[1]*data$prec + pars[2]*data$temp + pars[3]*data$mat_sum
  LU = pars[8]*data$total_fires * exp(pars[4]*data$ts_fire) + pars[5]*data$pasture + pars[6]*data$other_annual + pars[7]*data$other_perennial
  k = E + LU
  # Prediction of the model
  return ( 1-(exp(-(k))))
}

G <- function(pars) {
  # Extract parameters of the model
  LU = pars[1] * data$total_fires + pars[2] * data$ts_fire + pars[3] * data$pasture + 
  pars[4] * data$other_perennial + pars[5] * data$other_annual 
  E = pars[1]*data$prec + pars[2]*data$temp 
  k = E + LU
  # Prediction of the model
  return ( 1-(exp(-(k)))) #normalized_df$Bmax *   ) 
}

G <- function(pars) {
  # Extract parameters of the model
  k = pars[1] * data$age
  # Prediction of the model
  return (1-(exp(-(k))))
}

#pars = c(0.05,0.05, 0.05, 0.05,0.05,0.005,0.0005, 0.005)
pars = c(0.05, 0.05)
Gpred <- G(pars)
Gpred

# asymptote is being repeated.
# have so that the different selections have variance.
# normalize the variables so they dont vary wildly. divide the parameter values.

####################
#finding the right parameters

NLL = function(pars) {
  # Values prediced by the model
  if(pars[length(pars)] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  # before returning likelihood value, check if it is NA or inf or something problematic
  Gpred = G(pars)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = data$agbd - Gpred, mean = 0, sd = pars[length(pars)], log = TRUE), na.rm = TRUE)
  return(fun)
}

o = optim(par = pars, fn = NLL, hessian = FALSE)
print(o)

# meth0 = c("Nelder-Mead", "BFGS", "CG")
# for (i in meth0){
# o = optim(par = pars, fn = NLL, control = list(parscale = abs(pars)), 
#             hessian = FALSE, method = i)
# print(i)
# print(o)
# }

pred = G(o$par[1:length(o$par)])

outcome <- data.frame(data$agbd, pred)
outcome <- round(outcome, 3)

# Group by Column_A and calculate the median of Column_B for each group
median_values <- outcome %>%
  group_by(pred) %>%
  summarize(median_agbd = median(data.agbd, na.rm = TRUE))

head(outcome)
  
plot(median_values$median_agbd, median_values$pred, abline(0,1))

plot(outcome$data.agbd, outcome$pred, abline(0,1))


all_positive <- median_values %>%
  filter(pred >= 0)

tst <- lm(all_positive$pred ~ all_positive$median_agbd)
# 0.4716 R squared
# 0.202 R squared without LULC