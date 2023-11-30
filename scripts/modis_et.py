import ee
from geetools import batch
import pandas as pd

# Authenticate to Earth Engine
ee.Initialize()

def import_modis():
    r"""weeeee
    """

    modis = ee.ImageCollection("MODIS/061/MOD16A2")
    print(modis)
    
def yeee():
    print('ohno')


if __name__ == '__main__':
  import_modis()
  #yee()