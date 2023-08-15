####################################################################
########## Predicting forest regrowth from Mapbiomas data ##########
# first applied to the region of Paragominas, in the Western Amazon.
# Ana Avila - August 2023
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################

library(terra) # handling spatial data
library(tidyverse)
setwd("/home/aavila/forest_regrowth")

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
#santoro_raster <- mask(santoro_raster, regrowth)

# getting percent of mature forest cover within x neighboring patches:
mature_mask <- rast('./mapbiomas/mature_masks/0000000000-0000095232_mature_mask.tif')
mature_mask <- terra::crop(mature_mask, regrowth)
mature_mask[mature_mask > 0] <- 1
mature_mask <- subst(mature_mask, NA, 0)

range <- 21
mature_sum <- focal(mature_mask, range, sum, na.rm = TRUE) # run focal on these areas
mature_sum <- mask(mature_sum, regrowth) # select only the sum of areas surrounding secondary patches.

temp_cropped <- crop(temp, all_data)
prec_cropped <- crop(prec, all_data)
total_temp <- app(temp, sum)
total_prec <- app(prec, sum)

### save all data as a unified dataframe for easy visualization and modelling

all_data_santoro <- c(santoro_raster, regrowth, total_prec, total_temp, fire, last_LU)
file_name <- paste0(tempfile(), "_.tif")
# the stack is too big to convert right away
lu_tile <- makeTiles(all_data, c(2000,2000), file_name, na.rm = TRUE, overwrite=TRUE) 
lu_tile <- lapply(lu_tile, rast)

rast_to_df <- function(raster){
  coords <- crds(raster, df=FALSE, na.rm=TRUE)
  values_stack <- terra::extract(raster, coords, cells=FALSE, method="simple")
  central_df <- values_stack[complete.cases(values_stack), ]
  return(central_df)
}

data <- lapply(lu_tile, rast_to_df)
data <- bind_rows(data)
colnames(data) <- c('agbd', 'age', 'prec', 'temp', 'total_fires', 'ts_fire', 'last_LU')
#saveRDS(data, 'santoro_ESA_alldata.rds')

###############################################################

data <- readRDS('santoro_ESA_alldata.rds')


data <- data[data$last_LU %in% c(15, 41, 48),]
data$last_LU <- factor(data$last_LU)
dummy_LU <- as.data.frame(model.matrix(~ data$last_LU - 1))
names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')
data <- data[,-7]

minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
#data <- as.data.frame(lapply(data, minMax))

data_scaled <- cbind(agbd=data$agbd, scale(data[,2:ncol(data)]))
data_maxmin <- cbind(agbd=data$agbd, minMax(data[,2:ncol(data)]))

#data <- cbind(data, dummy_LU)
data <- as.data.frame(data_scaled)

######################################################################
#################        passing into the model     ##################
######################################################################

# why is everything converging to the asymptote?
# capture the temporal trend. is the temporal trend different everywhere?
# is the fitted asymptote the mean? that's the null model.



G <- function(pars) {
  E = pars['temp'] * data$temp + pars['prec'] * data$prec
  LU = pars['total_fires'] * data$total_fires + pars['ts_fire'] * data$ts_fire + pars['pasture'] * data$pasture + 
      pars['other_perennial'] * data$other_perennial + pars['other_annual'] * data$other_annual 
  k = E
  pars['B_0'] + pars['A'] * (1 - exp(-k))
}

pars = c(B_0 = 10, A = 100, temp =- 0.002, prec = 0.000005, total_fires = 0.05, 
ts_fire = 0.05, pasture = 0.05, other_perennial = 0.05, other_annual = 0.05, sd = 0.05)
Gpred <- G(pars)
head(Gpred)

NLL = function(pars) {
if(pars['sd'] < 0){ #avoiding NAs by keeping the st dev positive
  return(-Inf)
}
# Negative log-likelihood 
-sum(dnorm(x = data$agbd - Gpred, mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE)
}


o = optim(par = pars, fn = NLL, hessian = FALSE)   #, method = 'L-BFGS-B')
o

pred = G(o$par)

outcome <- data.frame(data$agbd, pred)
outcome <- round(outcome, 3)

median_values <- outcome %>%
  group_by(pred) %>%
  summarize(median_agbd = median(data.agbd, na.rm = TRUE))

plot(median_values$pred, median_values$median_agbd, abline(0,1), xlim=c(0, 100))

#####################################

fit <- lm(data$agbd ~ data$age)
summary(fit)
plot(data$agbd ~ data$age)
abline(fit)

sds <- aggregate(agbd ~ age, data, sd)
means <- aggregate(agbd ~ age, data, median)
sum_stats <- cbind(means, sds[,2])
colnames(sum_stats) <- c('age', 'median', 'sd')

ggplot(sum_stats,                               # ggplot2 plot with means & standard deviation
       aes(x = age,
           y = median)) + 
  geom_errorbar(aes(ymin = median - sd,
                    ymax = median + sd)) +
  geom_point() + theme(text = element_text(size = 20))  
