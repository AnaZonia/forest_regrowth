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

def get_minmax(img, roi):
    r"""get min and max values of the first image in an ImageCollection"""
    min_max_values = img.reduceRegion(reducer=ee.Reducer.minMax(), geometry=roi.geometry(), maxPixels = 4e10)
    return("min_max_values", min_max_values.getInfo())

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
    r"""Get mean evapotranspiration for MODIS over the whole time period across the basin.
    """
    roi = import_main_data(roi = "br_amaz")["roi"]
    modis = (ee.ImageCollection("MODIS/061/MOD16A2GF").filterDate("2002-01-01", "2021-12-01").filterBounds(roi.geometry()).select('ET')) # MODIS yearly ET ('min_max_values', {'ET_max': 735, 'ET_min': 8})
    def clip_image(image):
        return image.clip(roi)
    modis = modis.map(clip_image)
    modis_mean = modis.mean()
    return(modis_mean)


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
# def lulc():
#   last observed land use type - given secondary ages, get    
# if age is 1 in 2020, it was not forest in 2019
# which means I want category in 2019
# if age is 35 in 2020, it was not forest in 1985
# which means I want category in 1985

# transform secondary forest to mask
# mask land use at the correct band index to the secondary forest mask to the correct band index

########## LAND USE ##########
# VARIABLES OBTAINED
# number of years under each land use tnrype
# time since last observation of each land use type

# INDEX ## 3 = forest
# 15 = pasture
# 39 = soy
# 46 = coffee
# 20 = sugar cane
# 41 = other annual crop
# 48 = other perennial crop

# total years under each land use type
# calc_total_yrs <- function(masked_brick, val){
#   masked_brick[masked_brick != val] <- 0
#   total_past_years <- sum(masked_brick)
#   return(total_past_years/val)
# }

# pasture_total <- calc_total_yrs(lulc_brick_masked, 15)
# soy_total <- calc_total_yrs(lulc_brick_masked, 39)
# coffee_total <- calc_total_yrs(lulc_brick_masked, 46)
# sugar_total <- calc_total_yrs(lulc_brick_masked, 20)
# other_annual_total <- calc_total_yrs(lulc_brick_masked, 41)
# other_perennial_total <- calc_total_yrs(lulc_brick_masked, 48)


# time since specific land use type
# last_LU = lulc.eq(39).reduce(ee.Reducer.lastNonNull())


def fire():
    r""" Intakes MAPBIOMAS burned area and MAPBIOMAS fire frequency.
    Outputs Number of fires, time since last burn, and years since each fire."""
    fire = import_main_data(roi = "br_amaz")["fire"]
    roi = import_main_data(roi = "br_amaz")["roi"]
    # select only 2019 to 1985 (as we only want forests that burned before 2020)
    fire = fire.select(fire.bandNames().removeAll(["burned_coverage_2020", "burned_coverage_2021", "burned_coverage_2022"])
                                    .reverse())

    # fire has the value of the land use type that burned.
    # Transforming into a fire mask:
    fire = fire.gt(0)
    #display(fire)

    num_fires = fire.reduce('sum')
    print(get_minmax(num_fires, roi))

# --------------------------------------------------------------------------------------------

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


if __name__ == "__main__":
    #get_minmax(img)
    import_main_data(roi = "br_amaz")
    # smoothen_edges()
    ####### Climatic data
    # import_modis()
    # mcwd()
    # get_yearly_cwd()
    ####### LULC and Fire
    fire()
    ####### Export data
    # export_image_by_ecoregion(eco_id)
