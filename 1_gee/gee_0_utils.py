"""
Google Earth Engine Exploration Base Functions

This script provides a set of base functions and classes for exploring and processing
data in Google Earth Engine. It includes functionality for authentication, project
configuration, image export, map visualization, and land use/land cover data processing.

Key components:
- Earth Engine initialization
- Project configuration management
- Image export utility
- Map visualization tools
- Land use/land cover data processing

Author: [Your Name]
Date: [Current Date]
"""

import ee
import geemap

def initialize():
    """Initialize and authenticate Earth Engine."""
    try:
        ee.Initialize()
    except Exception as e:
        ee.Authenticate()
        ee.Initialize()

class ProjectConfig:
    """Singleton class to manage project configuration."""
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialize()
        return cls._instance

    def _initialize(self):
        """Initialize project configuration parameters."""
        self.data_folder = "projects/amazon-forest-regrowth/assets"
        self.roi = ee.FeatureCollection(f"{self.data_folder}/raw/biomes_br").geometry().dissolve()
        self.first_year = 1985
        self.last_year = 2020
        self.range_1986_2019 = range(self.first_year + 1, self.last_year)
        self.range_1985_2019 = range(self.first_year, self.last_year)
        self.range_1985_2020 = range(self.first_year, self.last_year + 1)

def export_image(img, name, folder=None, scale=None, crsTransform=None):
    """Export an Earth Engine image to an asset."""
    config = ProjectConfig()
    path = f"{config.data_folder}/{folder}" if folder else config.data_folder

    task_args = {
        'image': img,
        'description': name,
        'assetId': f"{path}/{name}",
        'region': config.roi,
        'crs': "EPSG:4326",
        'maxPixels': 4e12,
    }

    if scale is not None:
        task_args['scale'] = scale
    if crsTransform is not None:
        task_args['crsTransform'] = crsTransform

    task = ee.batch.Export.image.toAsset(**task_args)
    task.start()

class MapManager:
    """Manage map visualization for Earth Engine data."""
    def __init__(self):
        """Initialize map visualization parameters."""
        self.nicfi = ee.ImageCollection('projects/planet-nicfi/assets/basemaps/americas')
        self.basemap = self.nicfi.filter(ee.Filter.date('2020-03-01', '2020-12-01')).first()
        self.vis_planet = {'bands': ['R', 'G', 'B'], 'min': 64, 'max': 5454, 'gamma': 1.8}
        self.agbd_palette = {'min': 0, 'max': 400, 'palette': ['yellow', 'green']}
        self.age_palette = {'min': 0, 'max': 35, 'palette': ['white', 'red']}
        self.mask_palette = {'min': 0, 'max': 1, 'palette': ['white', 'black']}

    def init_map(self):
        """Initialize a new map centered on a specific location."""
        return geemap.Map(center=[-3.5, -60], zoom=10)

    def add_basemap(self, map_obj):
        """Add a NICFI basemap to the map."""
        map_obj.addLayer(self.basemap, self.vis_planet, 'NICFI Basemap')

    def add_agbd_layer(self, map_obj, agbd_data):
        """Add an Above Ground Biomass Density layer to the map."""
        map_obj.addLayer(agbd_data, self.agbd_palette, 'AGBD')

    def add_age_layer(self, map_obj, age_data):
        """Add an age layer to the map."""
        map_obj.addLayer(age_data, self.age_palette, 'Age')

def remap_band(band_name, img):
    """Remap a band to desired land use categories."""
    desired_values = ee.List([3, 6, 15, 39, 20, 40, 62, 41, 46, 47, 35, 48, 9, 21])
    mask_all_ones = ee.List.repeat(1, desired_values.size())
    band = img.select(ee.String(band_name))
    new_band = band.remap(desired_values, mask_all_ones, 0)
    return new_band.rename(ee.String(band_name))

def desired_lulc():
    """Process and return desired land use/land cover data."""
    config = ProjectConfig()

    age = ee.Image("projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_secondary_vegetation_age_v2").select("secondary_vegetation_age_2020")

    lulc = (ee.Image("projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_integration_v1")
            .select([f"classification_{year}" for year in config.range_1985_2020])
            .byte()
            .rename([str(year) for year in config.range_1985_2020])
            .updateMask(age))

    fire = (ee.Image("projects/mapbiomas-public/assets/brazil/fire/collection3/mapbiomas_fire_collection3_annual_burned_coverage_v1")
            .select([f"burned_coverage_{year}" for year in config.range_1985_2019])
            .byte()
            .rename([str(year) for year in config.range_1985_2019])
            .updateMask(age))

    remapped_image = lulc.bandNames().map(lambda band_name: remap_band(band_name, lulc))
    remapped_image = ee.ImageCollection(remapped_image).toBands()
    desired_mask = remapped_image.reduce("sum").eq(lulc.bandNames().size().getInfo())

    age = age.updateMask(desired_mask).rename("age")
    lulc = lulc.updateMask(desired_mask)
    fire = fire.updateMask(desired_mask)

    return age, lulc, fire
