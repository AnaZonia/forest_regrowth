
##################################################################################################

# The issue with matching UTM coordinates directly is a lot of points are lost. what we are plotting is actually a very small percentage of what is actually matching.
# with raster::resample, we can use nearest-neighbor.

GEDI = readRDS('GEDI_midamazon_dfunified.rds')
# GEDI = readRDS("0000000000-0000095232_GEDI.rds") for some reason, the method doesn't work with this file?
temp <- readRDS('./worldclim_dataframes/temp.rds')
temp <- temp[,c(1,2,37)]

# convert into raster
regrowth_raster <- rasterFromXYZ(regrowth_cleaned[,c(1,2,4)], crs = "+init=epsg:4326")
temp_raster <- rasterFromXYZ(temp[,c(1,2,3)], crs = "+init=epsg:4326")

# GEDI and soil data is irregular and can't be converted directly into a regular raster.

# making a raster from irregular data, using another raster as reference of size and resolution.
# df must be in form [lon;lat;data].
rasterFromXYX_irr <- function(df, ref_raster){   # dataframe to be converted, reference raster
  df = soil
  ref_raster = temp_raster
  e <- extent(min(df$lon), max(df$lon), min(df$lat), max(df$lat))
  # set up an 'empty' raster
  r <- raster(e, ncol=ncol(ref_raster), nrow=nrow(ref_raster))

  #GEDI <- GEDI[,c(1,3,2,4,5,6,7,8)]
  df[,3] <- as.factor(df[,3])
  x <- rasterize(df[, 1:2], r, df[,3])
  #x <- rasterize(GEDI[, 2:3], r, GEDI[,5])
  x

}












central_df <- cbind(agb_forest_age, last_burn = fire[match(agb_forest_age$xy,fire$xy),c("last_burn")])
central_df <- cbind(central_df, num_fires = fire[match(central_df$xy,fire$xy),c("num_fires")])
central_df <- cbind(central_df, pasture = lulc[match(central_df$xy,lulc$xy),c("pasture")])
central_df <- cbind(central_df, soy = lulc[match(central_df$xy,lulc$xy),c("soy")])
central_df <- cbind(central_df, other_perennial = lulc[match(central_df$xy,lulc$xy),c("other_perennial")])
central_df <- cbind(central_df, other_annual = lulc[match(central_df$xy,lulc$xy),c("other_annual")])
central_df <- cbind(central_df, tavg = temp[match(central_df$xy,temp$xy),c("mean")])
central_df <- cbind(central_df, prec = prec[match(central_df$xy,prec$xy),c("mean")])
central_df$last_burn[is.na(central_df$last_burn)] <- 1

central_df <- lapply(central_df, as.numeric)







################# passing into the model ##################

# pars[2]*tavg+pars[3]*prec+
# central_df$tavg, central_df$prec, 
#  0.001, 0.0001, 

columns = c(other_perennial, other_annual, pasture, soy, forest_age, num_fires, last_burn)

G = function(pars, columns) {
  # Extract parameters of the model
  Gmax = pars[1] #asymptote
  for (i in 1:length(columns)){
    k = pars[2]*other_perennial+pars[3]*pasture+pars[4]*soy+pars[5]*other_annual + pars[6]*other_annual + pars[7]*num_fires^(pars[8]*last_burn)
  }
  # Prediction of the model
  Gmax * (1 - exp(-k*(forest_age)))
}

####################
#finding the right parameters
par0 = c(50, 0.0001,  0.0001,  0.001, 0.001,  0.001, 0.1, 0.1)

Gpred = G(par0, central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

# dat$tavg, dat$prec, 

NLL = function(pars, dat) {
  # Values prediced by the model
  if(pars[9] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  Gpred = G(pars, dat$other_perennial, dat$other_annual, dat$pasture, dat$soy, dat$forest_age, dat$num_fires, dat$last_burn)
  #print(Gpred)
  # Negative log-likelihood 
  fun = -sum(dnorm(x = dat$agbd, mean = Gpred, sd = pars[9], log = TRUE))
  #print(pars)
  return(fun)
}

par0 = c(50, 0.001, 0.0001,  0.0001,  0.0001,  0.0001,  0.0001, 0.1, 0.1)
NLL(par0, dat=central_df)

o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
           hessian = FALSE, method = "CG")
print(o)

meth0 = c("Nelder-Mead", "BFGS", "CG")
for (i in meth0){
  o = optim(par = par0, fn = NLL, dat = central_df, control = list(parscale = abs(par0)), 
             hessian = FALSE, method = i)
  print(i)
  print(o)
}

#  central_df$tavg, central_df$prec, 

pred = G(o$par[1:9], central_df$other_perennial, central_df$other_annual, central_df$pasture, central_df$soy, central_df$forest_age, central_df$num_fires, central_df$last_burn)

plot(central_df$agbd, pred, abline(0,1))

