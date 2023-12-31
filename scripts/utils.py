import ee
import geemap

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
'''

def import_data(region):
    if region == "br_amazon":
        roi = ee.FeatureCollection("projects/ee-ana-zonia/assets/br_biomes").filter(ee.Filter.eq("id", 18413)).geometry()
    elif region == "br":
        roi = ee.FeatureCollection("projects/ee-ana-zonia/assets/br_shapefile").geometry()
    else:
        lulc = ee.Image("projects/mapbiomas-raisg/public/collection1/mapbiomas_raisg_panamazonia_collection1_integration_v1").byte()
    
    # Load images from MapBiomas Collection 8 for Land Use Land Cover and Burned Area
    if not region == "panamaz":
      lulc = (ee.Image("projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_integration_v1")
          .clip(roi) # clip data to region of interest specified above
          .select([f"classification_{year}" for year in range(first_year, last_year+1)])).byte()
      fire = (ee.Image("projects/mapbiomas-workspace/public/collection7_1/mapbiomas-fire-collection2-annual-burned-coverage-1")
          .clip(roi)
          .select([f"burned_coverage_{year}" for year in range(first_year, last_year)])).byte()


def export_grid_cell(feature):
    id = feature['id']
    grid_cell = fishnet.filter(ee.Filter.eq('id', id))
    img_export = img.clip(grid_cell)
    
    task = ee.batch.Export.image.toDrive(
        image = img_export,
        description = 'age_agbd',
        folder = 'fishnet_tiles',
        scale = 30,
        region = grid_cell.geometry(),
        maxPixels = 1e11
    )
    task.start()

def export_image_as_grid(img):
    fishnet = geemap.fishnet(roi, h_interval=2.0, v_interval=2.0, delta=0.5)
    # Set 'id' as a property of the features
    fishnet = fishnet.map(lambda feature: feature.set('id', feature.id()))
    fishnet_info = fishnet.getInfo()['features']
    img = img
    
    fishnet_info.map(export_grid_cell)