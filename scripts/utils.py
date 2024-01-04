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

def export_image(img, name, scale):
    # Create the export task
    task = ee.batch.Export.image.toAsset(
        image = img,
        description = f'{img}',
        assetId = f'projects/ee-ana-zonia/assets/{img}',
        scale = scale,
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
    

