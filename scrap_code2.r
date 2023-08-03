
### mature forest confusions

# mean biomass of surrounding mature forests around 3km (~100 pixels of 30m)
# biomass_range <- 101 # window size
# mature_biomass <- mask(santoro_raster, mature_mask)  # get only biomass of mature forests
# mature_total_biomass <- focal(mature_biomass, biomass_range, sum, na.rm = TRUE)

# Load required libraries
library(terra)
library(parallel)

tiles <- makeTiles(trim(mature_biomass), c(200,200), na.rm = TRUE, overwrite=TRUE)  # Adjust nx and ny to desired block size
tiles <- lapply(tiles, rast)

# function intakes buffer_size in meters and mean_biomass_in_x_pixels is self explanatory
calc_buffer_and_sum <- function(tile, buffer_size, mean_biomass_in_x_pixels) {

  mature_total_biomass <- focal(mature_biomass_msk, mean_biomass_in_x_pixels, mean, na.rm = TRUE)
  mature_biomass_msk <- mask(mature_total_biomass, rast(tile))

}

buffer_mask <- as.mask(dist)
mask_raster <- rast(r, mask = dist)

library(parallel)
tiles_processed <- mclapply(tiles, focal, 51, mean, na.rm = TRUE, ncores = 20)
saveRDS(tiles_processed, file = ".rds")


# Step 5: Combine the results from tiles into the final raster
result <- aggregate(tiles_processed, fun = sum)

plot(trim(mature_biomass))






data <- readRDS('santoro_ESA_alldata.rds')

data <- data[data$last_LU %in% c(15, 41, 48),]
data$last_LU <- factor(data$last_LU)
dummy_LU <- as.data.frame(model.matrix(~ data$last_LU - 1))
names(dummy_LU) <- c('pasture', 'other_annual', 'other_perennial')
data <- data[,-7]

minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
data <- as.data.frame(lapply(data, minMax))

#data <- as.data.frame(scale(data))

data <- cbind(data, dummy_LU)

######################################################################
#################        passing into the model     ##################
######################################################################

G <- function(pars) {
  E = pars['temp'] * data$temp + pars['prec'] * data$prec
  LU = pars['total_fires'] * data$total_fires + pars['ts_fire'] * data$ts_fire + pars['pasture'] * data$pasture + 
      pars['other_perennial'] * data$other_perennial + pars['other_annual'] * data$other_annual 
  k = E
  (1 - exp(-k))
}

pars = c(B_0 = 10, A = 100, temp = 0.5, prec = 0.5, total_fires = 0.05, ts_fire = 0.05, pasture = 0.05, other_perennial = 0.05, other_annual = 0.05,  sd = 0.05)
Gpred <- G(pars)
Gpred

  NLL = function(pars) {
  if(pars['sd'] < 0){ #avoiding NAs by keeping the st dev positive
    return(-Inf)
  }
  Gpred = G(pars)
  # Negative log-likelihood 
  -sum(dnorm(x = data$agbd - Gpred, mean = 0, sd = pars['sd'], log = TRUE), na.rm = TRUE)
  }


o = optim(par = pars, fn = NLL, hessian = FALSE)
o

pred = G(o$par[1:length(o$par)])

outcome <- data.frame(data$agbd, pred)
outcome <- round(outcome, 3)

median_values <- outcome %>%
  group_by(pred) %>%
  summarize(median_agbd = median(data.agbd, na.rm = TRUE))

plot(median_values$pred, median_values$median_agbd, abline(0,1))
