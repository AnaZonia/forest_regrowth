library(terra)
library(gstat)
library(sf)
library(jsonlite)


mat <- read.csv("data/mature_amazon_500m.csv")
coordinates(sec) <- ~ longitude + latitude
sec <- read.csv("data/unified_data.csv")

shp <- vect("data/br_shapefile")
shp <- st_as_sf(shp)
grid <- shp %>%
    st_make_grid(cellsize = 1, what = "centers") %>% # grid of points
    st_intersection(shp) # only within the polygon
grid <- as_Spatial(grid)
grid <- SpatialPixelsDataFrame(grid, data = data.frame(id = 1:length(grid)))
class(grid)

grid <- data.frame(grid)
names(grid) <- c("longitude", "latitude")


coords_df <- data.frame()
# Parse the JSON strings and extract the coordinates
for (x in mat$.geo) {
    coords <- fromJSON(x)$coordinates
    coords_df <- rbind(coords_df, c(longitude = coords[1], latitude = coords[2]))
}
names(coords_df) <- c("longitude", "latitude")
# Combine the original dataframe with the coordinates dataframe
mat <- cbind(mat, coords_df)


plot(var_mat <- variogram(b1 ~ 1, data = mat, cutoff = 7, width = 0.5), pl = TRUE)

plot(var_mat <- variogram(b1 ~ 1, data = mat, cutoff = 7, width = 0.25), pl = TRUE)

mod <- vgm(psill = 2700, model = "Exp", range = 2.25, nugget = 500)
plot(var_mat, mod2)

mod2 <- fit.variogram(var_mat, mod)

OK_est <- krige(b1 ~ 1, locations = ~longitude + latitude, data = mat, newdata = grid, model = mod2)
summary(OK_est)
pred.plot <- ggplot(aes(x = longitude, y = latitude), data = OK_est)
pred.plot <- pred.plot + geom_tile(aes(fill = var1.pred)) # change to var1.var for variance
pred.plot <- pred.plot + scale_fill_gradient(low = "yellow", high = "blue")
pred.plot + coord_equal()

OK_est <- krige(b1 ~ 1, locations = ~ longitude + latitude, data = mat, newdata = grid, model = mod2)

universal_est <- krige(b1 ~ longitude + latitude, locations = ~ longitude + latitude, data = mat, newdata = grid, model = mod2)

# edge effects with kriging - anywhere near the borders or the ocean
# run it by ecoregion

# start with low resolution grid

# universal kriging since there is a north-south trend in the data.
# block kriging for things that are not points

# co-krig

# maybe normalize the data 0-1

# cross-validate\

cross <- krige.cv(b1 ~ 1, locations = ~ longitude + latitude, data = mat, model = mod2, nfold = nrow(mat))

x <- 0
j <- 1
for (i in 1:25) {
  x <- x + ((j+2)^2 - j^2) * (25/i)^2
  j <- j + 2
}
x

