set.seed(123)

setwd("/home/aavila/forest_regrowth")

all_data_csv <- readRDS('santoro_ESA_alldata.rds')

head(all_data_csv)

x <- all_data_csv$age
y <- all_data_csv$agbd

G <- function(pars) {
  pars['a'] * (1 - exp(-pars['b'] * x))
}

neg_log_likelihood <- function(pars) {
  -sum(dnorm(y, G(pars), sd = 1, log = TRUE))
}

o = optim(par = c(a = 1, b = 1), fn = neg_log_likelihood)
o

pred = G(o$par[1:length(o$par)])

outcome <- data.frame(y, pred)
outcome <- round(outcome, 3)


##########################


# Set the seed for reproducibility
set.seed(123)

# Generate fake data
x <- seq(0, 10, by = 0.1)
y <- 100 * (1 - exp(-0.2 * x))^3 + rnorm(length(x), sd = 5)

# Define the Chapman-Richards growth model
chapman_richards <- function(x, a, b, c) {
  a * (1 - exp(-b * x))^c
}

# Define the negative log-likelihood function
NLL <- function(pars) {
  # Extract the parameters
  a <- pars[1]
  b <- pars[2]
  c <- pars[3]
  sd <- pars[4]

  # Compute the model predictions
  y_pred <- chapman_richards(x, a, b, c)

  # Compute the negative log-likelihood
  -sum(dnorm(y, mean = y_pred, sd = sd, log = TRUE))
}

# Set the initial values for the parameters
pars_init <- c(a = 10, b = 0.01, c = 1, sd = 5)

# Maximize the negative log-likelihood using nlm
fit <- nlm(NLL, p = pars_init)

# Print the estimated parameters
fit$estimate





# Load required libraries
library(terra)
library(parallel)

tiles <- makeTiles(trim(regrowth), c(200,200), na.rm = TRUE, overwrite=TRUE)  # Adjust nx and ny to desired block size
tiles <- lapply(tiles, rast)
length(tiles)
# function intakes buffer_size in meters and mean_biomass_in_x_pixels is self explanatory
calc_buffer_and_sum <- function(tile, buffer_size, mean_biomass_in_x_pixels) {
  buffer_size = 300
  r <- tiles[[1000]]
  pnt = as.points(r, values = FALSE)
  pol = buffer(pnt, 10)
  r[dist <= buffer_size] <- 1
  plot(dist)
  mature_biomass_crop <- crop(mature_biomass, dist)
  mature_biomass_msk <- mask(mature_biomass_crop, dist)
  mature_total_biomass <- focal(mature_biomass_msk, mean_biomass_in_x_pixels, mean, na.rm = TRUE)
  mature_biomass_msk <- mask(mature_total_biomass, rast(tile))

}

buffer_mask <- as.mask(dist)
mask_raster <- rast(r, mask = dist)


tiles_processed <- mclapply(tiles, calc_buffer_and_sum, 1000, 51, ncores = 10)

# Step 5: Combine the results from tiles into the final raster
result <- aggregate(tiles_processed, fun = sum)

plot(trim(mature_biomass))


