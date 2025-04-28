
# Tutorial at https://ecmwf-projects.github.io/copernicus-training-c3s/projections-cmip6.html

# CMIP6 climate projections https://cds.climate.copernicus.eu/datasets/projections-cmip6?tab=download
# Agroclimatic indicators derived https://cds.climate.copernicus.eu/datasets/sis-agroclimatic-indicators?tab=download


# General libs for file paths, data extraction, etc
from glob import glob
from pathlib import Path
from os.path import basename
import zipfile # To extract zipfiles
import urllib3 
urllib3.disable_warnings() # Disable warnings for data download via API
import os
import cftime
from xarray.coders import CFDatetimeCoder

# CDS API
import cdsapi

# Libraries for working with multi-dimensional arrays
import numpy as np
import xarray as xr
import pandas as pd

# Libraries for plotting and visualising data
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature

client = cdsapi.Client()

DATADIR = "./0_data/CMIP6"

experiments = ['historical', 'ssp126', 'ssp245', 'ssp585']

models = ['hadgem3_gc31_ll', 'inm_cm5_0', 'inm_cm4_8', 'ipsl_cm6a_lr', 
          'miroc_es2l', 'mpi_esm1_2_lr', 'ukesm1_0_ll']#, 'miroc6', 'canesm5']

variables = {
    # "ta": "air_temperature"
    # "pr": "precipitation",
    # "rsds": "surface_downwelling_shortwave_radiation",
    # "mrsos": "moisture_in_upper_portion_of_soil_column"
    "tas": "near_surface_air_temperature"
    # "huss": "near_surface_specific_humidity"
}


def download_data(experiment, start_year, end_year, models, variable="air_temperature"):
    """Download CMIP6 data from CDS for given experiment, years, models and variable."""
    for model in models:
        for start in range(start_year, end_year, 10):
            end = min(start + 11, end_year)  # Ensure no year exceeds end_year
            year_range = [str(y) for y in range(start, end)]

            try:
                # Attempt to download data
                client.retrieve("projections-cmip6", {
                    "temporal_resolution": "monthly",
                    "experiment": experiment,
                    "variable": variable,
                    "level": ["1000"],
                    "model": model,
                    "year": year_range,
                    "month": [f"{m:02d}" for m in range(1, 13)],  # All months
                    "area": [5.3, -74, -35, -34]
                }).download(f'{DATADIR}/{variable}/{variable}_{experiment}_{start}-{end-1}_{model}.zip')
                
                print(f"✅ Downloaded {variable} for {model} ({experiment}, {start}-{end-1})")

            except Exception as e:
                print(f"⚠️ Failed to download {variable} for {model} ({experiment}, {start}-{end-1}): {e}")



