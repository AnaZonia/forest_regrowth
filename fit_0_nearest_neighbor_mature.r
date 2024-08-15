################################################################################
#                                                                              #
#                 Get Nearest Mature forest                                           #
#                                                                              #
################################################################################

library(terra)
library(doParallel)

set.seed(0)

# Set up parallel processing
registerDoParallel(cores = 15)

#-------- SWITCHES

use_dist <- TRUE

#-------- FUNCTIONS

get_coords <- function(raster) {
    matrix <- as.matrix(raster, wide = TRUE)
    # matrix[matrix == 0] <- NA
    non_na_indices <- which(!is.na(matrix), arr.ind = TRUE)

    # Sample 10 random coordinates from non-NA values
    non_na_indices <- non_na_indices[sample(nrow(non_na_indices), 1000000, replace = FALSE), ]
    # should take about 6h to run 2.5 million rows with 20 cores

    # Extract sampled distances
    values <- matrix[cbind(non_na_indices[, 1], non_na_indices[, 2])]

    # Convert matrix indices to coordinates
    xy_coords <- xyFromCell(raster, cellFromRowCol(raster, non_na_indices[, 1], non_na_indices[, 2]))

    # Create data frame
    sampled_points_df <- data.frame(
        longitude = xy_coords[, 1],
        latitude = xy_coords[, 2],
        x_index = non_na_indices[, 2], #column number
        y_index = non_na_indices[, 1], #row number
        distance = values
    )
    return(sampled_points_df)
}

mean_biomass_within_buffer <- function(sec_pixel) {

    x <- sec_pixel[["x_index"]]
    y <- sec_pixel[["y_index"]]
    distance <- sec_pixel[["distance"]]

    # Calculate buffer limits
    buffer_radius <- ceiling(distance / 500) + 5

    # Define the extents of the buffer zone
    x_min <- max(1, x - buffer_radius)
    x_max <- min(ncol(mature_biomass), x + buffer_radius)
    y_min <- max(1, y - buffer_radius)
    y_max <- min(nrow(mature_biomass), y + buffer_radius)

    # Extract the buffer zone from the matrix
    buffer_zone <- mature_biomass[y_min:y_max, x_min:x_max]

    # Calculate the mean biomass
    non_zero_values <- buffer_zone[buffer_zone != 0]
    mean_biomass <- mean(non_zero_values, na.rm = TRUE)
    mean_biomass
    return(mean_biomass)
}

## -------- MAIN SCRIPT

mature_biomass <- rast("./data/mature_biomass_500m_countrywide.tif")

dist <- rast("./data/distance_1000m_countrywide.tif")
data <- get_coords(dist)

start_time <- Sys.time()
print(start_time)

nearest_biomass <- foreach(i = 1:nrow(data), .combine = "c") %dopar% {
    if (i %% 10000 == 0) {
        print(paste("Processing iteration:", i))
        print(paste("Time so far: ", as.numeric(difftime(
            Sys.time(), start_time,
            units = "mins"
        )), " minutes"))
    }
    mean_biomass_within_buffer(data[i, ])
}

print("nearest_biomass finished")

final_csv <- cbind(data, nearest_biomass)
final_csv <- na.omit(final_csv)
final_csv <- final_csv %>% rename("mature_biomass" = "nearest_biomass")

write.csv(final_csv, "./data/dist_mature_1000m_countrywide.csv", row.names = FALSE)
print("csv exported")

mat_biomass_raster <- rast(final_csv[, c("longitude", "latitude", "distance", "nearest_biomass")], type = "xyz", crs = crs(mature_biomass))

writeRaster(mat_biomass_raster, "./data/nearest_mat_biomass_1000m_countrywide.tif")
print("raster exported")
