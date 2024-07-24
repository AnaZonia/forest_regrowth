library(terra)
library(doParallel)

set.seed(0)
mature_biomass <- rast("./data/mature_biomass_1k.tif")

#-------- SWITCHES
use_dist <- TRUE

#-------- FUNCTIONS
get_coords <- function(raster) {
    matrix <- as.matrix(raster, wide = TRUE)
    matrix[matrix == 0] <- NA
    non_na_indices <- which(!is.na(matrix), arr.ind = TRUE)

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
        values = values
    )
    return(sampled_points_df)
}

mean_biomass_within_buffer <- function(sec_pixel, mature) {

    x <- sec_pixel[["x_index"]]
    y <- sec_pixel[["y_index"]]
    distance <- sec_pixel[["values"]]

    # Calculate buffer limits
    buffer_radius <- ceiling(distance / 500)

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

## -------- MAIN SCRIPT

if (use_dist) {
    # datafiles <- list(
    #     "5y_LULC",
    #     "10y_LULC",
    #     "15y_LULC",
    #     "all_LULC"
    # )

    # for (file in datafiles){
    # file = datafiles[1]
    # data <- read.csv(paste0("./data/", file, "_dist_amaz_500.csv"))
    # sec_vect <- vect(data, geom = c("longitude", "latitude"), crs = crs(mature_biomass))

    # template_raster <- rast(extent = ext(mature_biomass), resolution = res(mature_biomass), crs = crs(mature_biomass))
    # distance_raster <- rasterize(sec_vect, template_raster, field = "distance")
    # distance_df <- get_coords(distance_raster)
    # }


    dist <- rast("./data/distance.tif")
    data <- get_coords(dist)

    # Set up parallel processing
    registerDoParallel(cores = 15)

    nearest_biomass <- foreach(i = 1:nrow(data), .combine = "c") %dopar% {
        print(i)
        mean_biomass_within_buffer(data[i, ], mature_biomass)
    }

    final_csv <- cbind(data, nearest_biomass)
    write.csv(final_csv, "dist_mature_500m.csv")

    # write.csv(final_csv, paste0("./data/", datafile, "_mature_asymptote.csv"))

    mat_biomass_raster <- rast(final_csv[, c(1, 2, 5, 6)], type = "xyz", crs = crs(sec_dist))

    writeRaster(mat_biomass_raster, "./data/nearest_mat_biomass_500m.tif")

}else{
    mature_df <- get_coords(mature_biomass)
    sec_vect <- vect(data, geom = c("longitude", "latitude"), crs = crs(mature_biomass))
    mature_vect <- vect(mature_df, geom = c("x", "y"), crs = crs(mature_biomass))
    nn_mat <- nearest(sec_vect, mature_vect)
    save(nn_mat, file = "nn_mat.Rdata")
}


