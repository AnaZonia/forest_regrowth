import ee
import geemap
from functools import partial

# Authenticate to Earth Engine
try:
  ee.Initialize()
except Exception as e:
  ee.Authenticate()
  ee.Initialize(project='ee-ana-zonia')

'''
Choose region and import corresponding shapefile. Can be:
        br_amazon -> the Amazon biome in Brazilian territory
        br -> entire Brazil
        panamaz -> entire Amazonia, including territory across all Amazonian countries.
        
        *** panamaz has less land use categories and data processing is still unfinished for this data ***
        "projects/mapbiomasraisg/public/collection1/mapbiomas_raisg_panamazonia_collection1_integration_v1"
'''

def export_image(img, name):
    # Create the export task
    task = ee.batch.Export.image.toAsset(
        image = img,
        description = f'{name}',
        assetId = f'projects/ee-ana-zonia/assets/{name}',
        scale = 30,
        crs = 'EPSG:4326',
        maxPixels = 4e12
    )
    # Start the export task
    task.start()

def map_image(img, min, max):
    vis = {
        'min': min,
        'max': max,
        'palette': ['blue', 'red'],
    }
    
    Map = geemap.Map()
    Map.addLayer(img, vis)
    return Map
    
# Deleting assets in bulk from projects. Taken from https://gis.stackexchange.com/questions/467363/batch-deleting-of-earth-engine-assets

asset_list = ee.data.listAssets("projects/ee-ana-zonia/assets")["assets"]
asset_list

def conditional_asset_rm(x, starts_with):
    """Deletes asset if starts with starts_with """
    id = x["id"]              # users/username/file  or projects/project-name/assets/file
    findex = 5 if id.startswith("users") else 3
    name = x["name"]          # projects/earthengine-legacy/assets/users/username/file or projects/project-name/assets/file
    f = name.split("/")[findex]  # file
    if (f.startswith(starts_with)):
        ee.data.deleteAsset(id)
        return f"Deleted asset {id}"
    
    # return 0


# [conditional_asset_rm(x, "unif") for x in asset_list]


