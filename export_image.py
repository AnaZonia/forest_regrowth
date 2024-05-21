def export_image(img, name):
    # Create the export task
    task = ee.batch.Export.image.toAsset(
        image = img,
        description = f'{name}',
        assetId = f'projects/amazon-forest-regrowth/assets/{name}',
        scale = 30,
        crs = 'EPSG:4326',
        maxPixels = 4e12
    )
    # Start the export task
    task.start()