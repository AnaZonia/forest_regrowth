import matplotlib.pyplot as plt

tst = age_agbd.select('age').updateMask(land_use_15_years.select("last_LU")).updateMask(one_hectare_mask)

# 19, 24, 29
histogram_ee = mature_biomass.reduceRegion(
    geometry = tst.geometry(), reducer = ee.Reducer.histogram(), maxPixels=1e13
)

hist = histogram_ee.getInfo()

# list(hist.keys())
# Get the bin centers and bin counts
bin_centers = hist['mature_biomass']['bucketMeans']
bin_counts = hist['mature_biomass']['histogram']

# Plot the histogram
plt.bar(bin_centers, bin_counts, width=1)
plt.xlabel('Age')
plt.ylabel('Count')
plt.title('Histogram of Ages with 5 years of land use history or less')
plt.show()
