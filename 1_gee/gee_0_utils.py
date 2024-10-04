import ee
import geemap

def initialize():
    # Authenticate to Earth Engine
    try:
        ee.Initialize()
    except Exception as e:
        ee.Authenticate()
        ee.Initialize()


class ProjectConfig:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialize()
        return cls._instance

    def _initialize(self):
        self.data_folder = "projects/amazon-forest-regrowth/assets"
        self.roi = ee.FeatureCollection(f"{self.data_folder}/raw/biomes_br").geometry().dissolve()
        self.first_year = 1985
        self.last_year = 2020
        self.range_1986_2019 = range(self.first_year + 1, self.last_year)
        self.range_1985_2019 = range(self.first_year, self.last_year)
        self.range_1985_2020 = range(self.first_year, self.last_year + 1)




def export_image(img, name, folder = None, scale = None, crsTransform = None):
    config = ProjectConfig()

    path = config.data_folder
    if folder:
        path = f"{path}/{folder}"

    # Create the export task
    task_args = {
        'image': img,
        'description': f"{name}",
        'assetId': f"{path}/{name}",
        'region': config.roi,
        'crs': "EPSG:4326",
        'maxPixels': 4e12,
    }

    # Add scale or crsTransform to the task arguments
    if scale is not None:
        task_args['scale'] = scale
    if crsTransform is not None:
        task_args['crsTransform'] = crsTransform

    # Create the export task
    task = ee.batch.Export.image.toAsset(**task_args)
    
    # Start the export task
    task.start()



class MapManager:
    def __init__(self):
        self.nicfi = ee.ImageCollection('projects/planet-nicfi/assets/basemaps/americas')
        self.basemap = self.nicfi.filter(ee.Filter.date('2020-03-01', '2020-12-01')).first()
        self.vis_planet = {'bands': ['R', 'G', 'B'], 'min': 64, 'max': 5454, 'gamma': 1.8}
        self.agbd_palette = {'min': 0, 'max': 400, 'palette': ['yellow', 'green']}
        self.age_palette = {'min': 0, 'max': 35, 'palette': ['white', 'red']}
        self.mask_palette = {'min': 0, 'max': 1, 'palette': ['white', 'black']}

    def init_map(self):
        return geemap.Map(center=[-3.5, -60], zoom=10)

    def add_basemap(self, map_obj):
        map_obj.addLayer(self.basemap, self.vis_planet, 'NICFI Basemap')

    def add_agbd_layer(self, map_obj, agbd_data):
        map_obj.addLayer(agbd_data, self.agbd_palette, 'AGBD')

    def add_age_layer(self, map_obj, age_data):
        map_obj.addLayer(age_data, self.age_palette, 'Age')



"""
## Removing pixels with undesired land use categories

Some land use categories are not relevant to the model (such as rocky surfaces or mangroves)

All pixels with **at least one observation of the undesired land use history** are used to make a mask, to leave behind only pixels with occurrences of only desired land use types.


Land use types we are interested in:

    3 = forest
    15 = pasture
    39 = soy
    20 = sugar cane
    21 = mosaic of uses
    40 = rice
    62 = cotton
    41 = other temporary crop
    46 = coffee
    47 = citrus
    35 = palm oil
    48 = other perennial crop
    9 = forest plantationantation
"""

# For each band, convert pixels with desired land use types to 1 - undesired types to zero
def remap_band(band_name, img):
    # List the categories that are DESIRED to be maintained
    desired_values = ee.List([3, 6, 15, 39, 20, 40, 62, 41, 46, 47, 35, 48, 9, 21])
    mask_all_ones = ee.List.repeat(1, desired_values.size())

    band = img.select(ee.String(band_name))
    new_band = band.remap(desired_values, mask_all_ones, 0)
    return new_band.rename(ee.String(band_name))

def desired_lulc():
    config = ProjectConfig()

    # import ages from MapBiomas
    age = ee.Image(
        "projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_secondary_vegetation_age_v2"
    ).select("secondary_vegetation_age_2020")

    # Load images from MapBiomas Collection 8 for Land Use Land Cover and Burned Area
    lulc = (
        ee.Image(
            "projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_integration_v1"
        )
        .select([f"classification_{year}" for year in config.range_1985_2020])
        .byte()
        .rename([str(year) for year in config.range_1985_2020])
    ).updateMask(age)

    fire = (
        ee.Image(
            "projects/mapbiomas-public/assets/brazil/fire/collection3/mapbiomas_fire_collection3_annual_burned_coverage_v1"
        )
        .select([f"burned_coverage_{year}" for year in config.range_1985_2019])
        .byte()
        .rename([str(year) for year in config.range_1985_2019])
        .updateMask(age)
    )


    # Map the function over the band names
    remapped_image = lulc.bandNames().map(lambda band_name: remap_band(band_name, lulc))
    # make mask by adding all pixels that add up to the total number of years (all pixels with desired categories)
    remapped_image = ee.ImageCollection(remapped_image).toBands()
    desired_mask = remapped_image.reduce("sum").eq(lulc.bandNames().size().getInfo())

    age = age.updateMask(desired_mask).rename("age")
    lulc = lulc.updateMask(desired_mask)
    fire = fire.updateMask(desired_mask)

    return age, lulc, fire

