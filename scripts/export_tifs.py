import ee
from geetools import batch

# Authenticate to Earth Engine
ee.Initialize()

# Load the images and feature collections
age = ee.Image('users/celsohlsj/public/secondary_vegetation_age_collection71_v5').select('classification_2020')
biomass = ee.Image('projects/ee-ana-zonia/assets/biomass_2020')
sd = ee.Image('projects/ee-ana-zonia/assets/biomass_sd_2020')
cwd = ee.Image('projects/ee-ana-zonia/assets/cwd_chave')

# Regions of interest
amazon_biome = ee.FeatureCollection('projects/ee-ana-zonia/assets/amazon_biome_border')
indig_land = ee.FeatureCollection('projects/ee-ana-zonia/assets/indig_land')
ecoregions = ee.FeatureCollection("RESOLVE/ECOREGIONS/2017")
ecoregions = ecoregions.filterBounds(amazon_biome.geometry())
ecoregions_list = ecoregions.toList(ecoregions.size())

# Clip to Amazon biome
age = age.clip(amazon_biome).updateMask(age.gt(0))
biomass = biomass.clip(amazon_biome)
sd = sd.clip(amazon_biome)
cwd = cwd.clip(amazon_biome)

# Reproject to 10m
biomass_10m = biomass.reproject(crs=biomass.projection(), scale=10)
cwd_10m = cwd.reproject(crs=age.projection(), scale=10)
sd_10m = sd.reproject(crs=age.projection(), scale=10)

# Reaggregate to 30m (mean value)
aggregated_biomass = biomass_10m.reduceResolution(reducer=ee.Reducer.mean()).reproject(crs=age.projection())
aggregated_cwd = biomass_10m.reduceResolution(reducer=ee.Reducer.mean()).reproject(crs=age.projection())
aggregated_sd = biomass_10m.reduceResolution(reducer=ee.Reducer.mean()).reproject(crs=age.projection())

# Mask only to regions with age greater than zero (secondary forests)
aggregated_biomass = aggregated_biomass.updateMask(age).float()
aggregated_sd = aggregated_sd.updateMask(age).float()
aggregated_cwd = aggregated_cwd.updateMask(age).float()

def export_data(img, prefix, ecoregion):
    # Export the image to an Earth Engine asset.
    task = ee.batch.Export.image.toDrive(
        image=img,
        description=f"{prefix}_{ecoregion.get('ECO_ID').getInfo()}",
        folder='drive_export',
        region=ecoregion.geometry(),
        crs=img.projection().getInfo()['crs'],
        crsTransform=img.projection().getInfo()['transform'],
        skipEmptyTiles = TRUE,
        maxPixels=4e10
    )
    task.start()

for i in range(31):
    ecoreg = ee.Feature(ecoregions_list.get(i))
    img_clipped = age.clip(ecoreg)
    cwd_clipped = aggregated_cwd.clip(ecoreg)
    sd_clipped = aggregated_sd.clip(ecoreg)
    
    export_data(cwd_clipped, 'aggregated_cwd', ecoreg)
    export_data(img_clipped, 'age', ecoreg)
    export_data(sd_clipped, 'aggregated_sd', ecoreg)

