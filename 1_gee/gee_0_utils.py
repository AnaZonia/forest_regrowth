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

Author: Ana Catarina Avila
Date: 2024-08-29
"""

import ee
import geemap

# ------------------------------ Initialization ------------------------------

def initialize():
    """Initialize and authenticate Earth Engine."""
    try:
        ee.Initialize()
    except Exception as e:
        ee.Authenticate()
        ee.Initialize()

# ------------------------------ Project Configuration ------------------------------

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

# ------------------------------ Image Export ------------------------------

def export_image(img, name, region, folder = None, scale = None, crsTransform = None):
    """Export an Earth Engine image to an asset."""
    config = ProjectConfig()
    path = f"{config.data_folder}/{folder}" if folder else config.data_folder

    task_args = {
        'image': img,
        'description': name,
        'assetId': f"{path}/{name}",
        'region': region,
        'crs': "EPSG:4326",
        'maxPixels': 4e12,
    }

    if scale is not None:
        task_args['scale'] = scale
    elif crsTransform is not None:
        task_args['crsTransform'] = crsTransform

    task = ee.batch.Export.image.toAsset(**task_args)
    task.start()




# ------------------------------ Land Use/Land Cover Data Processing ------------------------------

# The MapBiomas Collection 9 land use/land cover data is mapped to the following classes:

# 3 <- Forest Formation
# 6 <- Floodable Forest 
# 15 <- Pasture
# 20 <- Sugar Cane
# 21 <- Mosaic of Uses
# 35 <- Palm Oil
# 39 <- Soybean
# 40 <- Rice
# 41 <- Other temporary crops
# 46 <- Coffee
# 47 <- Citrus
# 48 <- Other perennial crops
# 62 <- Cotton

# 9 <- Forest Plantation (excluded)

# Note: These classes have different accuracies, and the accuracy is expected to decrease further in the past.

def remap_band(band_name, img):
    """Remap a band to desired land use categories."""
    desired_values = ee.List([3, 6, 15, 20, 21, 35, 39, 40, 41, 46, 47, 48, 62])
    mask_all_ones = ee.List.repeat(1, desired_values.size())
    band = img.select(ee.String(band_name))
    new_band = band.remap(desired_values, mask_all_ones, 0)
    return new_band.rename(ee.String(band_name))

def desired_lulc():
    """Process and return desired land use/land cover data."""
    config = ProjectConfig()

    age = ee.Image("projects/mapbiomas-public/assets/brazil/lulc/collection9/mapbiomas_collection90_secondary_vegetation_age_v1").select("secondary_vegetation_age_2020")

    lulc = (ee.Image("projects/mapbiomas-public/assets/brazil/lulc/collection9/mapbiomas_collection90_integration_v1")
            .select([f"classification_{year}" for year in config.range_1985_2020])
            .byte()
            .rename([str(year) for year in config.range_1985_2020])
            .updateMask(age))

    remapped_image = lulc.bandNames().map(lambda band_name: remap_band(band_name, lulc))
    remapped_image = ee.ImageCollection(remapped_image).toBands()
    # select only pixels with exclusively the desired land use histories (exclude all instances of land use types we are not interested in)
    desired_mask = remapped_image.reduce("sum").eq(lulc.bandNames().size().getInfo())

    age = age.updateMask(desired_mask).rename("age")
    lulc = lulc.updateMask(desired_mask)

    return age, lulc

# ------------------------------ Create grid ------------------------------
# Used in gee_5_mature and gee_6_write_csv

def create_grid(image, region_name, cell_size = 10000, file_name = None):
    
    config = ProjectConfig()
    biomes = ee.Image(f"{config.data_folder}/categorical").select("biome")

    pixels_to_sample = biomes.updateMask(image)

    if region_name == "amazon":
        biome_index = 1
    elif region_name == "atlantic":
        biome_index = 4

    roi = ee.FeatureCollection(f"{config.data_folder}/raw/biomes_br").filterMetadata('CD_Bioma', 'equals', biome_index).geometry()

    # First, sample locations based only on the age band
    grid = geemap.create_grid(roi, cell_size, 'EPSG:4326')

    # Function to sample one point per valid cell
    def sample_cell(cell):
        sampled_fc = pixels_to_sample.stratifiedSample(
            numPoints = 1,
            classBand = 'biome',
            region = cell.geometry(),
            scale = cell_size,
            geometries = True,
            dropNulls = True
        )

        # Only return a feature if we found one
        return ee.Feature(ee.Algorithms.If(
            sampled_fc.size().gt(0),
            sampled_fc.first(),
            # Return a placeholder that we can filter out later
            ee.Feature(ee.Geometry.Point([0, 0])).set('is_null', True)
        ))

    samples = grid.map(sample_cell)

    # Filter out placeholder features before exporting
    samples = samples.filter(ee.Filter.notEquals('is_null', True))

    if file_name is None:
        return samples
    else:
        export_name = f"grid_{cell_size//1000}k_{region_name}_{file_name}"

        export_task = ee.batch.Export.table.toAsset(
            collection = samples,
            description = export_name,
            assetId = f"{config.data_folder}/{export_name}"
        )
        export_task.start()
