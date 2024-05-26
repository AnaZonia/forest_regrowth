# def prec_anomaly(month_prec):
#     anom = month_prec.subtract(mean_prec).divide(sd_prec)
#     drought_mask = anom.lt(0)
#     return anom.updateMask(drought_mask)

# monthly_anom = prec.map(prec_anomaly)

# def anom_yearly(year):
#     # get anomaly of the driest month of the year
#     anom_year = monthly_anom.filter(ee.Filter.calendarRange(year, year, 'year')).reduce(ee.Reducer.min())
#     return anom_year.set('year', year)

# yearly_anom = ee.ImageCollection.fromImages(years.map(anom_yearly)).toBands()

# new_band_names = ['yearly_anom_{}'.format(year) for year in yearlist]
# yearly_anom = yearly_anom.rename(new_band_names)


# anom = prec.first().subtract(mean_prec).divide(sd_prec)
# drought_mask = anom.lt(0)
# plot = anom.updateMask(drought_mask)


# years = ee.List.sequence(1985, 2019)

# terraclim = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE') \
#               .filterDate('1985-01-01', '2019-12-31') \
#               .map(lambda image : image.clip(roi)) \
#               .map(lambda image: image.reduceResolution(ee.Reducer.median(), bestEffort=True, maxPixels=1024) \
#                                        .reproject(crs='EPSG:4326', scale=10000))

# maxtemp = terraclim.select('tmmx').map(lambda image: image.multiply(0.1))
# mintemp = terraclim.select('tmmn').map(lambda image: image.multiply(0.1))
# radiation = terraclim.select('srad').map(lambda image: image.multiply(0.1))

# def mean_yearly(year):
#     mean_year = radiation.filter(ee.Filter.calendarRange(year, year, 'year')).reduce(ee.Reducer.mean())
#     return mean_year.set('year', year)

# mean_rad = ee.ImageCollection.fromImages(years.map(mean_yearly)).toBands()

# yearlist = range(1985, 2020) # Generate a list of years from 1985 to 2019
# new_band_names = ['rad_{}'.format(year) for year in yearlist] # Append 'si_' to each year
# mean_rad = mean_rad.rename(new_band_names)

# vis = {
#     'min': 0,
#     'max': 5000,
#     'palette': ['blue', 'red'],
# }

# Map = geemap.Map(center=[-10, -40], zoom=4)
# Map.addLayer(mean_rad.select('rad_2010'), vis, 'rad')
# Map


# MODIS Evapotranspiration

# Calculate mean AET from MODIS and calculate CWD from the precipitation values as in Celso

# Since MODIS only goes back to 2000, for now we are stuck with a fixed value for ET

# Function to calculate mean AET and add year property
# select for mature forests since the values can be put off by deforestation (causes lower ET)

start = "2002-01-01"
end = "2019-12-31"

modis = (
    ee.ImageCollection("MODIS/061/MOD16A2GF")
    .filterDate(start, end)
    .select("ET", "ET_QC")
    .map(lambda image: image.clip(roi))
    .map(
        lambda image: image.reduceResolution(
            ee.Reducer.median(), bestEffort=True, maxPixels=1024
        ).reproject(crs="EPSG:4326", scale=10000)
    )
)


# code sourced from https://spatialthoughts.com/2021/08/19/qa-bands-bitmasks-gee/
def bitwise_extract(input, from_bit, to_bit):
    mask_size = ee.Number(1).add(to_bit).subtract(from_bit)
    mask = ee.Number(1).leftShift(mask_size).subtract(1)
    return input.rightShift(from_bit).bitwiseAnd(mask)


def apply_QA_mask(image):
    QA = image.select("ET_QC")
    ET = image.select("ET").multiply(
        0.0125
    )  # multiply by the scale 0.1, divide by 8 to get daily info
    cloud_mask = bitwise_extract(QA, 3, 4).lte(0)
    qa_mask = bitwise_extract(QA, 5, 7).lte(1)
    mask = cloud_mask.And(qa_mask)
    return ET.updateMask(mask).set("system:time_start", image.get("system:time_start"))


# mask quality of pixels
modis_masked = modis.map(apply_QA_mask).map(lambda image: image.updateMask(mature_mask))


# Loop through the months and filter images
def monthly_et(month):
    eight_day_images = modis_masked.select("ET").filter(
        ee.Filter.calendarRange(month, month, "month")
    )
    month_et = (
        eight_day_images.median()
        .multiply(30)
        .reduceRegions(ecoregions, reducer=ee.Reducer.median(), scale=10000)
    )
    month_et = month_et.reduceToImage(["month_et"], ee.Reducer.first())
    return month_et.set("month", month)


monthly_et_imgcol = ee.ImageCollection.fromImages(months.map(monthly_et))

yearlist = range(2002, 2019)  # Generate a list of years from 1985 to 2019
monthlist = range(1, 12)
new_band_names = []
for year in yearlist:
    for month in monthlist:
        new_band_names = [f"mcwd_{month}_{year}"]
monthly_et_img = monthly_et_imgcol.toBands().rename(new_band_names)
