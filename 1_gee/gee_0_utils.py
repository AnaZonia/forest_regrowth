import ee
import geemap

def initialize():
    # Authenticate to Earth Engine
    try:
        ee.Initialize()
    except Exception as e:
        ee.Authenticate()
        ee.Initialize()


def init_map():
    Map = geemap.Map()
    Map = geemap.Map(center=[-3.5, -60], zoom = 10)
    return Map


# Function to get data folder and ROI
def global_variables():
    data_folder = "projects/amazon-forest-regrowth/assets"
    roi = ee.FeatureCollection(f"{data_folder}/raw/biomes_br").geometry().dissolve()
    first_year = 1985
    last_year = 2020
    # 1986 - 2019, years included in analysis
    years = range((first_year + 1), last_year)
    return data_folder, roi, years, first_year, last_year


def export_image(img, name, scale):
    # Create the export task
    task = ee.batch.Export.image.toAsset(
        image=img,
        description=f"{name}",
        assetId=f"projects/amazon-forest-regrowth/assets/{name}",
        scale=scale,
        region = roi,
        crs="EPSG:4326",
        maxPixels=4e12,
    )
    # Start the export task
    task.start()
