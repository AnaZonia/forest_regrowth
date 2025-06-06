


library(terra)
library(utils)
library(dplyr)


# Settings
experiments <- c("historical", "ssp126", "ssp245", "ssp585") 

models <- c("hadgem3_gc31_ll",
    "inm_cm5_0", "inm_cm4_8",
    "ipsl_cm6a_lr", "miroc_es2l", "mpi_esm1_2_lr"#, "ukesm1_0_ll" has two types of files - clarify later if worth including
)

cmip6_dir <- "./0_data/CMIP6/"

variables <- list(
    "ta" = "air_temperature",
    "pr" = "precipitation",
    "rsds" = "surface_downwelling_shortwave_radiation",
    "mrsos" = "moisture_in_upper_portion_of_soil_column",
    "huss" = "near_surface_specific_humidity",
    "tas" = "near_surface_air_temperature"
)

# For precipitation conversion
seconds_in_month <- c(
    2678400, 2419200, 2678400, 2592000, 2678400, 2592000, 2678400,
    2678400, 2592000, 2678400, 2592000, 2678400
) # typical non-leap year

# For solar radiation conversion)
hours_in_month <- c(744, 672, 744, 720, 744, 720, 744, 744, 720, 744, 720, 744)




# ====================================
# FUNCTION: process_model_experiment
# Description: Processes each model × experiment × variable .nc file
# into a yearly-averaged raster stack
# ====================================

process_model_experiment <- function(experiment, model, var_dir, var_short) {

    model_dir <- file.path(var_dir, model)

    # List all relevant .nc files
    nc_files <- list.files(
        path = model_dir,
        pattern = paste0(".*", experiment, ".*\\.nc$"),
        full.names = TRUE
    )

    if (length(nc_files) == 0) {
        message("No .nc files found for experiment: ", experiment)
        return(NULL)
    }

    combined_raster <- NULL

    for (nc_file in nc_files) {
        r <- rast(nc_file)
        message("Loaded: ", nc_file)

        num_layers <- nlyr(r)
        num_chunks <- ceiling(num_layers / 12)
        yearly_rasters <- list()

        for (i in seq_len(num_chunks)) {
            start_layer <- ((i - 1) * 12) + 1
            end_layer <- min(i * 12, num_layers)
            chunk <- r[[start_layer:end_layer]]
            
            if (var_short %in% c("tas", "ta", "huss", "mrsos")) {
                # Variables that should be averaged
                layer <- app(chunk, fun = mean, na.rm = TRUE)
            } else if (var_short == "pr") {
                # Precipitation: weight by seconds and sum
                weighted <- chunk * seconds_in_month
                layer <- app(weighted, fun = sum, na.rm = TRUE)
            } else if (var_short == "rsds") {
                # Solar radiation: convert W/m² to kWh/m²/year
                # Each monthly layer (W/m²) * hours in month / 1000 = kWh/m²/month
                weighted_solar <- list()
                for (month_idx in 1:min(12, nlyr(chunk))) {
                    monthly_layer <- chunk[[month_idx]]
                    # Convert W/m² to kWh/m²: multiply by hours, divide by 1000
                    monthly_energy <- monthly_layer * hours_in_month[month_idx] / 1000
                    weighted_solar[[month_idx]] <- monthly_energy
                }
                # Sum all monthly energy values for annual total
                solar_stack <- do.call(c, weighted_solar)
                layer <- app(solar_stack, fun = sum, na.rm = TRUE)
            }
            yearly_rasters[[i]] <- layer
        }

        yearly_stack <- do.call(c, yearly_rasters)

        if (is.null(combined_raster)) {
            combined_raster <- yearly_stack
        } else {
            combined_raster <- c(combined_raster, yearly_stack)
        }
    }

    # Name bands based on the year range
    if (grepl("historical", basename(nc_files[[1]]))) {
        names(combined_raster) <- seq(1950, 2014)
    } else {
        names(combined_raster) <- seq(2015, 2074)
    }

    # Output file path
    output_file <- file.path(model_dir, paste0(var_short, "_", experiment, "_", model, ".tif"))

    writeRaster(
        combined_raster,
        filename = output_file,
        filetype = "GTiff",
        overwrite = TRUE
    )
    message("Saved: ", output_file)
}

# ====================================
# FUNCTION: standardize_extents
# Description: Aligns extent, resolution, and CRS of a list of rasters
# ====================================
standardize_extents <- function(experiment_rasters) {
    extents <- lapply(experiment_rasters, ext)

    xmin_ref <- max(unlist(lapply(extents, function(x) x$xmin)))
    xmax_ref <- min(unlist(lapply(extents, function(x) x$xmax)))
    ymin_ref <- max(unlist(lapply(extents, function(x) x$ymin)))
    ymax_ref <- min(unlist(lapply(extents, function(x) x$ymax)))

    min_extent <- ext(xmin_ref, xmax_ref, ymin_ref, ymax_ref)
    ref_raster <- experiment_rasters[[1]]

    # Reproject and crop all rasters to the smallest extent
    experiment_rasters_aligned <- lapply(experiment_rasters, function(r) {
        r <- resample(r, ref_raster) # resample to a common resolution
        r <- project(r, crs(ref_raster)) # reproject to the same CRS
        r <- crop(r, min_extent) # crop to the smallest extent
        return(r)
    })

    return(experiment_rasters_aligned)
}



# ====================================
# MAIN LOOP: Process each variable
# ====================================
for (var_short in names(variables)) {

    var_full <- variables[[var_short]]
    var_dir <- file.path(cmip6_dir, var_full)

    # Create one raster per model per experiment
    for (model in models) {
        # Unzip model-specific files
        zip_files <- list.files(
            path = var_dir,
            pattern = paste0(".*", model, ".*\\.zip$"),
            full.names = TRUE
        )
        
        if (length(zip_files) == 0) {
            message("No zip files found for model: ", model)
            next
        }

        model_dir <- file.path(var_dir, model)
        dir.create(model_dir, showWarnings = FALSE)

        for (zip_file in zip_files) {
            unzip(zipfile = zip_file, exdir = model_dir)
            message("Unzipped: ", zip_file)
        }

        # Process each experiment
        for (experiment in experiments) {
            process_model_experiment(experiment, model, var_dir, var_short)
        }
    }

    # ====================================
    # Combine rasters from all models for each experiment
    # ====================================
    for (experiment in experiments) {
        experiment_list <- c()

        for (model in models) {
            # Combine all models for each experiment
            model_files <- list.files(
                path = file.path(var_dir, model),
                pattern = paste0(var_short, "_", experiment, "_.*\\.tif$"),
                full.names = TRUE
            )
            experiment_list <- c(experiment_list, model_files)
        }

        if (length(experiment_list) == 0) {
            message("No model files found for experiment: ", experiment)
            next
        }

        experiment_rasters <- lapply(experiment_list, rast)

        experiment_rasters_aligned <- standardize_extents(experiment_rasters)

        s <- sds(experiment_rasters_aligned)
        combined_raster <- app(s, mean)

        # Set names based on the year range
        if (grepl("historical", experiment)) {
            names(combined_raster) <- seq(1950, 2014)
        } else {
            names(combined_raster) <- seq(2015, 2074)
        }

        output_file <- paste0(cmip6_dir, var_full, "_", experiment, ".tif")

        writeRaster(
            combined_raster,
            filename = output_file,
            filetype = "GTiff",
            overwrite = TRUE
        )
        message("Saved: ", output_file)
    }
}



