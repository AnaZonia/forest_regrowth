# Load necessary libraries
library(tidyverse)

# Load data
data <- read.csv("your_file.csv")  # replace "your_file.csv" with your actual file path

# Define grid and parameters
lambda <- 0.1  # sand settling rate, as specified
kappa <- 0.01  # diffusion coefficient, adjust as needed
dx <- 1  # grid spacing in x-direction
dy <- 1  # grid spacing in y-direction
dt <- 1  # time step for the simulation

# Initialize sand grid based on initial dumping site
# Assuming the CSV has columns: x, y, current_magnitude, current_direction, sand_initial
grid <- data %>%
  mutate(vx = current_magnitude * cos(current_direction),
         vy = current_magnitude * sin(current_direction),
         sand = sand_initial)  # Initialize sand amount at each grid cell

# Function to calculate advection term
advection <- function(sand, vx, vy, dx, dy, i, j, grid) {
  adv_x <- -vx * (grid$sand[i+1, j] - grid$sand[i-1, j]) / (2 * dx)
  adv_y <- -vy * (grid$sand[i, j+1] - grid$sand[i, j-1]) / (2 * dy)
  return(adv_x + adv_y)
}

# Function to calculate diffusion term
diffusion <- function(sand, kappa, dx, dy, i, j, grid) {
  diff_x <- (grid$sand[i+1, j] - 2 * grid$sand[i, j] + grid$sand[i-1, j]) / (dx^2)
  diff_y <- (grid$sand[i, j+1] - 2 * grid$sand[i, j] + grid$sand[i, j-1]) / (dy^2)
  return(kappa * (diff_x + diff_y))
}

# Function to update sand concentration at each time step
update_sand <- function(grid, lambda, kappa, dx, dy, dt) {
  updated_grid <- grid
  for (i in 2:(nrow(grid) - 1)) {
    for (j in 2:(ncol(grid) - 1)) {
      # Calculate advection and diffusion terms
      adv <- advection(grid$sand, grid$vx[i, j], grid$vy[i, j], dx, dy, i, j, grid)
      diff <- diffusion(grid$sand, kappa, dx, dy, i, j, grid)
      # Update sand amount with full model equation
      updated_grid$sand[i, j] <- grid$sand[i, j] +
        dt * (adv + diff - lambda * grid$sand[i, j])
    }
  }
  return(updated_grid)
}

# Time loop for diffusion model
timesteps <- 100  # Adjust based on required simulation time
for (t in 1:timesteps) {
  grid <- update_sand(grid, lambda, kappa, dx, dy, dt)
}

# View final sand distribution
print(grid)
