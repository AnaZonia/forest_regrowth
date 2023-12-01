#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generating data for "Predicting Forest Regrowth in Amazonia given Land Use History"
Author: Ana Catarina Avila
Date: Nov 2023
"""

import ee
import geemap
from geetools import batch
import pandas as pd

# Authenticate to Earth Engine
ee.Initialize()
# Initialize map
Map = geemap.Map(center=[-10, -40], zoom=4)
Map


def get_minmax(img, geom):
    r"""get min and max values of the first image in an ImageCollection"""
    min_max_values = img.first().reduceRegion(reducer=ee.Reducer.minMax(), geometry=geom)
    print("min_max_values", min_max_values.getInfo())

def import_main_data(roi):
    r""" Import images from ee.
    Clips to region of interest (roi):
        - roi = "br_amaz": The Brazilian Amazon
        - roi = "br": All of Brazil
        - roi = "panamaz": The entire Amazonia (under development)
    - removes urban areas
    Returns a dictionary containing:
        - FeatureCollections: roi, indig_land, ecoregions, soil
        - Images: age, biomass, lulc, fire
    """
    
    if roi == "br_amaz":
        roi = ee.FeatureCollection("projects/ee-ana-zonia/assets/amazon_biome_border").filter(ee.Filter.eq("id", 18413))
    elif roi == "br":
        roi = ee.FeatureCollection("projects/ee-ana-zonia/assets/br_shapefile")
    else:
        # land use for the pan amazonian dataset (spanning all amazonian countries)
        panamaz = ee.Image("projects/mapbiomas-raisg/public/collection1/mapbiomas_raisg_panamazonia_collection1_integration_v1")
        # to apply here method of Silva Junior et al 2020 to get secondary forest ages
        # note: there's less land use categories here than for the Brazilian territory.
    
    Map.addLayer(roi, {}, 'roi')
    # FeatureCollections
    indig_land = ee.FeatureCollection("projects/ee-ana-zonia/assets/indig_land").filterBounds(roi.geometry())
    ecoregions = ee.FeatureCollection("RESOLVE/ECOREGIONS/2017").filterBounds(roi.geometry())
    soil = soil = ee.FeatureCollection('projects/ee-ana-zonia/assets/DSMW').filterBounds(roi.geometry())
    
    # Images
    urban = ee.Image("DLR/WSF/WSF2015/v1")
    biomass = ee.Image("projects/ee-ana-zonia/assets/biomass_2020").clip(roi)
    lulc = ee.Image("projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_integration_v1").clip(roi)
    fire = ee.Image("projects/mapbiomas-workspace/public/collection7_1/mapbiomas-fire-collection2-annual-burned-coverage-1").clip(roi)

    # Load the images and feature collections
    if not roi == "panamaz":
        age = (ee.Image("users/celsohlsj/public/secondary_vegetation_age_collection71_v5")
            .select("classification_2020") # since the biomass data is from 2020
            .clip(roi))
        # include only pixels with age greater than zero (secondary forests)
        age = age.updateMask(age.gt(0))

    return {"roi": roi, "indig_land": indig_land, "ecoregions": ecoregions,
            "age": age, "soil": soil, "biomass":biomass, "lulc":lulc, "fire": fire}

def smoothen_edges():
    r"""
    This section is done to account for the "edge" pixels. To attribute a more accurate biomass value
    for a 30x30 age pixel that is at the edge of a 100x100 biomass pixel, I downsample biomass,
    average the values, and then reaggregate the 10x10 pixels to 30x30, realigning with age.
    """
    
    # add and rename bands
    img_export = age.addBands(biomass).addBands(sd).addBands(cwd).rename(["age", "agbd", "agbd_sd", "cwd"])
    # Reproject to 10m
    img_export = img_export.reproject(crs=age.projection(), scale=10)
    # Reaggregate to 30m (mean value)
    img_export = img_export.reduceResolution(reducer=ee.Reducer.mean()).reproject(crs=age.projection())
    # Mask only to regions with age greater than zero (secondary forests)
    img_export = img_export.updateMask(age).float()

    for eco_id in ecoregions.aggregate_array("ECO_ID").getInfo():
        export_image_by_ecoregion(eco_id)


"""
~~~~~~   Climatic Data   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def import_modis():
    r"""weeeee
    """
    roi = import_main_data(roi = "br_amaz")["roi"]
    modis = (ee.ImageCollection("MODIS/061/MOD16A2GF").select("ET")
             .filterBounds(roi).filterDate("2001-01-01", "2022-12-01")) # MODIS yearly ET

# def mcwd():

def get_yearly_cwd():
    # Referenced from https://github.com/celsohlsj/RasterMCWD/blob/v1.0.8/MCWD.R
    
    # Load TerraClimate dataset and filter by date
    terraclimate = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE") \
        .filter(ee.Filter.date("1985-07-01", "2020-08-01"))
    # Select the chosen band and clip to the Amazon biome - note there is scaling for some values.
    prec = terraclimate.select("pr").filterBounds(roi) # precipitation, 0 - 617 mm
    tmin = terraclimate.select("tmmn").filterBounds(roi) # min temp, 10.9 - 24 deg C
    tmax = terraclimate.select("tmmx").filterBounds(roi) # max temp, 21.3 - 36.4 deg C
    
    # get wd as rainfall - evapotranspiration
    # result = 

"""
~~~~~~   LULC and Fire   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""




# def fire():
""" """
# num_fires <- sum(fire_brick_masked)

# # find when was last fire observed (how many years before observed regrowth)
# fire_brick_flipped <- fire_brick_masked [[ c(rev(order(names(fire_brick_masked)))) ]]
# last_fire <- which.lyr(fire_brick_flipped == 1)


# def mature_layer():
""" """
# import mapbiomas lulc script
# select for the pixels between 200 and 300
# mature[mature > 300 | mature < 200] <- NA # only leave behind mature values
# mature_mask <- mask(mature, app(mature, fun = sum))
# consider only pixels that have been mature the whole time

# get biomass of surrounding mature forests
# with adjustable buffer size



"""
~~~~~~   Export Data   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

def export_image_by_ecoregion(eco_id):
    ecoreg =  ecoregions.filter(ee.Filter.eq("ECO_ID",eco_id))
    img = img_export.clip(ecoreg)
    proj = img.projection().getInfo()
    task = ee.batch.Export.image.toDrive(
            image=img,
            description=f"img_export_{eco_id}",
            folder="drive_export",
            region=ecoreg.geometry(),
            crs=proj["crs"],
            crsTransform=proj["transform"],
            skipEmptyTiles = True,
            maxPixels=4e10
        )
    task.start()

def ee_array_to_df(img, list_of_bands):
    r"""Transforms client-side ee.Image.getRegion array to pandas.DataFrame."""
    # Attempting to extract more than 1.048.576 values will result in an error.
    
    # arr = make array with getRegion from Image object
    
    df = pd.DataFrame(arr)

    # Rearrange the header.
    headers = df.iloc[0]
    df = pd.DataFrame(df.values[1:], columns=headers)

    # Remove rows without data inside.
    df = df[["longitude", "latitude", "time", *list_of_bands]].dropna()

    # Convert the data to numeric values.
    for band in list_of_bands:
        df[band] = pd.to_numeric(df[band], errors="coerce")

    # Convert the time field into a datetime.
    df["datetime"] = pd.to_datetime(df["time"], unit="ms")

    # Keep the columns of interest.
    df = df[["time","datetime",  *list_of_bands]]

    return df



if __name__ == "__main__":
    # get_minmax(img, geom)
    import_main_data(roi = "br_amaz")
    # smoothen_edges()
    ####### Climatic data
    import_modis()
    # mcwd()
    # get_yearly_cwd()
    ####### LULC and Fire
    ####### Export data
    # export_image_by_ecoregion(eco_id)
    # ee_array_to_df(img, list_of_bands)