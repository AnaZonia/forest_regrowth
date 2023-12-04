

# making masks
# 1.1 Anthropic, urban and Water Mask
empty = ee.Image().byte()

for i in range(first_year, last_year + 1):
    year = 'classification_' + str(i)
    anthropic = lulc.select(year).remap(
        [15, 19, 39, 20, 40, 62, 41, 36, 46, 47, 48, 9, 21, 24, 30],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        0
    ).rename(year)
    empty = empty.addBands(anthropic)

anthropic_mask = empty.select(empty.bandNames().slice(1)) # we only want areas that DO show anthropic activity this year

# Replace 'YOUR_WATER_IMAGE_ID' with the actual water image ID you are working with
w_mask  = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("max_extent").clip(roi).remap([0,1],[1,0]);

urban = ee.Image("DLR/WSF/WSF2015/v1").clip(roi)
inverse_mask = urban.eq(1).Not()  # This will invert the mask, city pixels will be False (or 0) and non-urban pixels will be True (or 1)
urban = urban.updateMask(inverse_mask)  # This will replace 1 values (city pixels) with NA
urban_mask = urban.unmask(1)  # This will replace NA values (non-urban pixels) with 1

# Define a color palette
# palette = ['blue', 'red']  # Blue for 0, Red for 1
# Add the layer to the map with the color palette
# Map.addLayer(w_mask, {'min': 0, 'max': 1, 'palette': palette}, 'w_mask')
# Map.addLayer(urban_mask, {'palette': palette}, 'urban_mask')
# Map.addLayer(anthropic_mask.select('classification_2017'), {'min': 0, 'max': 1, 'palette': palette}, "anthropic_mask")
# Map


# 1. Reclassifying MapBiomas Data # Step 1
empty = ee.Image().byte();

for i in range(first_year, last_year+1):
    year = 'classification_' + str(i)
    forest = lulc.select(year)
    forest = forest.remap([3,6], [1,1], 0) # Forest Formation and Flooded Forest classes from MapBiomas Project
    empty = empty.addBands(forest.rename(ee.String(year)))

mapbiomas_forest = empty.select(empty.bandNames().slice(1))


# 2. Mapping the Annual Increment of Secondary Forests # Step 2
regro = ee.Image().byte()
defor = ee.Image().byte()

for i in range(first_year, last_year):  # 1986-2020
    year1 = f'classification_{i}'
    year2 = f'classification_{i + 1}'
    a_mask = anthropic_mask.select(year1);
    forest1 = mapbiomas_forest.select(year1).remap([0, 1], [0, 2])  # Change forest pixels in 1985 to 2 years old
    forest2 = mapbiomas_forest.select(year2)
    # addition is 0 if was nonforest before and after; 1 if it was gained; 2 if it was forest before and then was lost; 3 if it was forest in both.
    sforest = forest1.add(forest2).multiply(a_mask).multiply(w_mask).multiply(urban_mask)
    for_gain = sforest.remap([0, 1, 2, 3], [0, 1, 0, 0]).rename(year2)
    for_loss = sforest.remap([0, 1, 2, 3], [0, 0, 1, 0]).rename(year2)
    regro = regro.addBands(for_gain)
    defor = defor.addBands(for_loss)

regro = regro.select(regro.bandNames().slice(1))  # Shows all years in which forest was gained.
# here, we could just mask by pixels that are forest in 2020 and find the year of last gain.

# 3. Mapping the Annual Extent of Secondary Forests # Step 3
extent = ee.Image().byte()
# add pixels that gained forest in 1986
extent = extent.addBands(regro.select('classification_1986').rename('classification_1986'))

for i in range(first_year + 1, last_year): #1987 to 2020
    year = f'classification_{i}' #1986
    year2 = f'classification_{i + 1}' #1987
    for_gain = regro.select(year2)
    acm_forest = extent.select(year).add(for_gain) #pixels that gained forest in 1986 + pixels that gained forest in 1987
    old_values = list(range(37))
    new_values = [0, 1] + [1] * 35
    remap = acm_forest.remap(old_values, new_values)
    # mask (multiply) by pixels that were shown to be forest in 1987, hence eliminating any that may have regrown in 1986 but lost cover in 1987
    extent = extent.addBands(remap.multiply(mapbiomas_forest.select(year2)).rename(year2))

extent = extent.select(extent.bandNames().slice(1))


# 4. Calculating and Mapping the Age of Secondary Forests # Step 4
ages = ee.Image().byte()
ages = ages.addBands(extent.select('classification_1986').rename('classification_1986'))
ages = ages.slice(1)


for i in range(first_year, last_year):  # 1986-2020
    year = f'classification_{i + 1}'
    sforest = sforest_ext.select(year)
    age_forest = empty.add(sforest)
    f_year = mapbiomas_forest.select(year)
    age_forest = age_forest.multiply(f_year)
    sforest_age = sforest_age.addBands(age_forest.rename(year))
    empty = age_forest

sforest_age = sforest_age.updateMask(sforest_age)
vizpar = {'min': 1, 'max': last_year - first_year, 'palette': ['blue', 'red']}  # Blue for 0, Red for 1
Map = geemap.Map(center=[-10, -40], zoom=4)
Map.addLayer(sforest_age.select('classification_1987'), vizpar, "classification_2017")
Map










#######################


# 3.1 Secondary Forest Loss
empty = ee.Image().byte()
empty2 = ee.Image().byte()
ext = sforest_all.select('classification_1986')
ext = ext.rename('classification_1986')
empty = empty.addBands(ext)

for i in range(first_year + 1, last_year):  # 1987-2020
    year = f'classification_{i}'
    year2 = f'classification_{i + 1}'
    sforest = sforest_all.select(year2)
    acm_forest = empty.select(year).add(sforest)
    old_values = list(range(39))
    new_values = [0, 1] + [1] * 37
    remap = acm_forest.remap(old_values, new_values)
    mask = mapbiomas_forest.select(year2).remap([0, 1], [500, 1])
    loss = remap.add(mask).remap([1, 2, 500, 501], [0, 0, 0, 1])
    empty2 = empty2.addBands(loss.rename(year2))
    empty = empty.addBands(remap.multiply(mapbiomas_forest.select(year2)).rename(year2))

sforest_loss = empty2.select(empty2.bandNames().slice(1))
# generates a mask with yearly secondary forest cover loss