library(terra)
library(doParallel)

set.seed(0)
mature_biomass <- rast("./data/mature_biomass_1k.tif")
dist <- rast("./data/distance.tif")

get_coords <- function(raster){
    matrix <- as.matrix(raster, wide = TRUE)
    matrix[matrix == 0] <- NA
    non_na_indices <- which(!is.na(matrix), arr.ind = TRUE)

    # # Sample 10 random coordinates from non-NA values
    # non_na_indices <- non_na_indices[sample(nrow(non_na_indices), 1000000, replace = FALSE), ]

    # Extract sampled distances
    values <- matrix[cbind(non_na_indices[, 1], non_na_indices[, 2])]

    # Convert matrix indices to coordinates
    xy_coords <- xyFromCell(raster, cellFromRowCol(raster, non_na_indices[, 1], non_na_indices[, 2]))

    # Create data frame
    sampled_points_df <- data.frame(
        x = xy_coords[, 1],
        y = xy_coords[, 2],
        xcoord = non_na_indices[, 2],
        ycoord = non_na_indices[, 1],
        value = sampled_distances
    )
    return(sampled_points_df)
}

mature_df <- get_coords(mature_biomass)

data <- read.csv("./data/all_LULC.csv")
sec_dist <- vect(data, geom = c("longitude", "latitude"), crs = crs(mature_biomass))
mature_vect <- vect(mature_df, geom = c("x", "y"), crs = crs(mature_biomass))

nn_mat <- nearest(sec_dist, mature_vect)

mean_biomass_within_buffer <- function(sec_pixel, mature) {
    x <- sec_pixel[[3]]
    y <- sec_pixel[[4]]
    dist <- sec_pixel[[5]]

    # Calculate buffer limits
    buffer_radius <- ceiling(dist / 500)
    
    # Define the extents of the buffer zone
    x_min <- max(1, x - buffer_radius)
    x_max <- min(ncol(mature), x + buffer_radius)
    y_min <- max(1, y - buffer_radius)
    y_max <- min(nrow(mature), y + buffer_radius)

    # Extract the buffer zone from the matrix
    buffer_zone <- mature[y_min:y_max, x_min:x_max]

    # Calculate the mean biomass
    non_zero_values <- buffer_zone[buffer_zone != 0]
    mean_biomass <- mean(non_zero_values, na.rm = TRUE)

    return(mean_biomass)
}

# Set up parallel processing
registerDoParallel(cores = 15)

nearest_biomass <- foreach(i = 1:nrow(sampled_points_df), .combine = "c") %dopar% {
    print(i)
    mean_biomass_within_buffer(sampled_points_df[i, ], mature_biomass)
}

final_csv <- cbind(sampled_points_df, nearest_biomass)

write.csv(final_csv, "./data/nearest_biomass.csv")

mat_biomass_raster <- rast(final_csv[,c(1,2,5,6)], type = "xyz", crs = crs(sec_dist))

writeRaster(mat_biomass_raster, "./data/mat_biomass_raster.tif")