
import ee
import geemap
import geopandas as gpd
from shapely.geometry import shape

ee.Authenticate()
ee.Initialize()
Map = geemap.Map()

ages_2020 = ee.Image('users/celsohlsj/public/secondary_vegetation_age_collection71_v5').select('classification_2020')
biomass_2020 = ee.Image('projects/ee-ana-zonia/assets/biomass_2020')
sd = ee.Image('projects/ee-ana-zonia/assets/biomass_sd_2020')
cwd = ee.Image('projects/ee-ana-zonia/assets/cwd_brazilian_amazon_chave')
biome_shp = ee.FeatureCollection('projects/ee-ana-zonia/assets/amazon_biome_border')

# Clip, mask, and rename
ages_2020 = ages_2020.clip(biome_shp.geometry()).updateMask(ages_2020.gt(0)).rename(['age'])
Map.addLayer(ages_2020, {}, 'ages_2020')

# Resample to 10m
biomass_10m = biomass.reproject(crs=ages_2020.projection(), scale=10)
cwd_10m = cwd.reproject(crs=ages_2020.projection(), scale=10)
sd_10m = sd.reproject(crs=ages_2020.projection(), scale=10)

# Mask only to regions with age greater than zero
biomass_10m = biomass_10m.updateMask(ages_2020)
sd_10m = sd_10m.updateMask(ages_2020)
cwd_10m = cwd_10m.updateMask(ages_2020)

# Reaggregate to 30m (mean value)
aggregated_biomass = biomass_10m.reduceResolution(
    reducer=ee.Reducer.mean(),
    maxPixels=1024
).reproject(crs=ages_2020.projection().crs(), scale=30)

aggregated_cwd = cwd_10m.reduceResolution(
    reducer=ee.Reducer.mean(),
    maxPixels=1024
).reproject(crs=ages_2020.projection().crs(), scale=30)

aggregated_sd = sd_10m.reduceResolution(
    reducer=ee.Reducer.mean(),
    maxPixels=1024
).reproject(crs=ages_2020.projection().crs(), scale=30)

# Get the covering grid
grid = ages_2020.geometry().coveringGrid(
    projection=ages_2020.projection(),
    size=100000
)

# Print the covering grid (you may want to use grid.getInfo() for a more readable output)
print(grid)


