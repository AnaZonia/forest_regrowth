import ee
from geetools import batch
import pandas as pd

# Authenticate to Earth Engine
ee.Initialize()

# Regions of interest
amazon_biome = ee.FeatureCollection('projects/ee-ana-zonia/assets/amazon_biome_border')
indig_land = ee.FeatureCollection('projects/ee-ana-zonia/assets/indig_land')
ecoregions = ee.FeatureCollection("RESOLVE/ECOREGIONS/2017")
ecoregions = ecoregions.filterBounds(amazon_biome.geometry())
ecoregions_list = ecoregions.toList(ecoregions.size())

# Load the images and feature collections
age = (ee.Image('users/celsohlsj/public/secondary_vegetation_age_collection71_v5')
       .select('classification_2020')
       .clip(amazon_biome))
age = age.updateMask(age.gt(0)) # include only pixels with age greater than zero (secondary forests)
biomass = ee.Image('projects/ee-ana-zonia/assets/biomass_2020').clip(amazon_biome)
sd = ee.Image('projects/ee-ana-zonia/assets/biomass_sd_2020').clip(amazon_biome)
cwd = ee.Image('projects/ee-ana-zonia/assets/cwd_chave').clip(amazon_biome)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this section is done to account for the "edge" pixels. To attribute a more accurate biomass value
# for a 30x30 age pixel that is at the edge of a 100x100 biomass pixel, I downsample biomass,
# average the values, and then reaggregate the 10x10 pixels to 30x30, realigning with age.

# add and rename bands
img_export = age.addBands(biomass).addBands(sd).addBands(cwd).rename(['age', 'agbd', 'agbd_sd', 'cwd'])
# Reproject to 10m
img_export = img_export.reproject(crs=age.projection(), scale=10)
# Reaggregate to 30m (mean value)
img_export = img_export.reduceResolution(reducer=ee.Reducer.mean()).reproject(crs=age.projection())
# Mask only to regions with age greater than zero (secondary forests)
img_export = img_export.updateMask(age).float()

def export_by_ecoregion(eco_id):
    ecoreg =  ecoregions.filter(ee.Filter.eq("ECO_ID",eco_id))
    img = img_export.clip(ecoreg)
    proj = img.projection().getInfo()
    task = ee.batch.Export.image.toDrive(
            image=img,
            description=f"img_export_{eco_id}",
            folder='drive_export',
            region=ecoreg.geometry(),
            crs=proj["crs"],
            crsTransform=proj["transform"],
            skipEmptyTiles = True,
            maxPixels=4e10
        )
    task.start()

for eco_id in ecoregions.aggregate_array("ECO_ID").getInfo():
    export_by_ecoregion(eco_id)


def ee_array_to_df(arr, list_of_bands):
    """Transforms client-side ee.Image.getRegion array to pandas.DataFrame."""
    df = pd.DataFrame(arr)

    # Rearrange the header.
    headers = df.iloc[0]
    df = pd.DataFrame(df.values[1:], columns=headers)

    # Remove rows without data inside.
    df = df[['longitude', 'latitude', 'time', *list_of_bands]].dropna()

    # Convert the data to numeric values.
    for band in list_of_bands:
        df[band] = pd.to_numeric(df[band], errors='coerce')

    # Convert the time field into a datetime.
    df['datetime'] = pd.to_datetime(df['time'], unit='ms')

    # Keep the columns of interest.
    df = df[['time','datetime',  *list_of_bands]]

    return df