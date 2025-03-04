
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


DATADIR = "./0_data/CMIP6"

experiments = ['historical', 'ssp126', 'ssp245', 'ssp585']

models = ['hadgem3_gc31_ll', 'inm_cm5_0', 'inm_cm4_8', 'ipsl_cm6a_lr', 
          'miroc_es2l', 'mpi_esm1_2_lr', 'ukesm1_0_ll', 'miroc6', 'canesm5']

variables = {
    "pr": "precipitation",
    "rsds": "surface_downwelling_shortwave_radiation",
    "mrsos": "moisture_in_upper_portion_of_soil_column",
    "tas": "near_surface_air_temperature"
    # "huss": "near_surface_specific_humidity"
}

client = cdsapi.Client()

def download_data(experiment, start_year, end_year, models, variable="air_temperature"):
    for j in models:
        for start in range(start_year, end_year, 10):
            end = min(start + 11, end_year)  # Ensure no year exceeds end_year
            year_range = [str(y) for y in range(start, end)]

            try:
                # Attempt to download data
                client.retrieve("projections-cmip6", {
                    "temporal_resolution": "monthly",
                    "experiment": experiment,
                    "variable": variable,
                    # "level": ["1000"],
                    "model": j,
                    "year": year_range,
                    "month": [f"{m:02d}" for m in range(1, 13)],  # All months
                    "area": [5.3, -74, -35, -34]
                }).download(f'{DATADIR}/{variable}/{variable}_{experiment}_{start}-{end-1}_{j}.zip')
                
                print(f"✅ Successfully downloaded {variable} for {j} ({experiment}, {start}-{end-1})")

            except Exception as e:
                print(f"⚠️ Skipping {variable} for {j} ({experiment}, {start}-{end-1}): {e}")



for variable_short, variable in variables.items():
    # DOWNLOAD DATA FOR HISTORICAL PERIOD
    download_data("historical", 1985, 2015, models, variable)

    # DOWNLOAD DATA FOR FUTURE SCENARIOS
    for i in experiments[1:]:
        download_data(i, 2015, 2050, models, variable)

    # Extract zip files

    cmip6_zip_paths = glob(f'{DATADIR}/{variable}/*.zip')
    for j in cmip6_zip_paths:
        with zipfile.ZipFile(j, 'r') as zip_ref:
            zip_ref.extractall(f'{DATADIR}/{variable}/')

    # Find all NetCDF files
    data_files = sorted(glob(f"{DATADIR}/{variable}/*.nc"))

    # Extract the base name (without year) of each file
    grouped_files = {}
    for file in data_files:
        base_name = "_".join(os.path.basename(file).split('_')[:-1])  # Remove year part
        grouped_files.setdefault(base_name, []).append(file)

    # Function to merge NetCDF files
    def merge_ncdf(files, output_path):
        datasets = [xr.open_dataset(f) for f in files]
        combined = xr.concat(datasets, dim="time")  # Merge along time dimension
        combined.to_netcdf(output_path)
        print(f"Saved merged dataset: {output_path}")

    # Merge files and save
    for base_name, files in grouped_files.items():
        output_file = f"{DATADIR}/{variable}/{base_name}_merged.nc"
        merge_ncdf(files, output_file)

    cmip6_nc = list()
    cmip6_nc_rel = glob(f'{DATADIR}/{variable}/*merged.nc')
    for i in cmip6_nc_rel:
        cmip6_nc.append(os.path.basename(i))


    # Load and prepare data
    # Function to aggregate in geographical lat lon dimensions

    def geog_agg(fn):
        ds = xr.open_dataset(f'{DATADIR}/{variable}/{fn}', engine='netcdf4')
        exp = ds.attrs['experiment_id']
        mod = ds.attrs['source_id']
        da = ds[f'{variable_short}']
        weights = np.cos(np.deg2rad(da.lat))
        weights.name = "weights"
        da_weighted = da.weighted(weights)
        da_agg = da_weighted.mean(['lat', 'lon'])
        da_yr = da_agg.groupby('time.year').mean()
        if variable in ['air_temperature', 'daily_maximum_near_surface_air_temperature']:
            da_yr = da_yr - 273.15
        da_yr = da_yr.assign_coords(model=mod)
        da_yr = da_yr.expand_dims('model')
        da_yr = da_yr.assign_coords(experiment=exp)
        da_yr = da_yr.expand_dims('experiment')
        da_yr.to_netcdf(path=f'{DATADIR}/{variable}/agg_{exp}_{mod}_{str(da_yr.year[0].values)}.nc')

    for i in cmip6_nc:
        try:
            geog_agg(i)
        except: print(f'{i} failed')

    data_ds = xr.open_mfdataset(f'{DATADIR}/{variable}/agg*.nc')
    data_ds.load()
    data = data_ds[f'{variable_short}']
    # data = data.squeeze("plev")

    data_90 = data.quantile(0.9, dim='model')
    data_10 = data.quantile(0.1, dim='model')
    data_50 = data.quantile(0.5, dim='model')

    fig, ax = plt.subplots(1, 1, figsize = (16, 8))

    colours = ['black','red','green','blue']
    for i in np.arange(len(experiments)):
        ax.plot(data_50.year, data_50[i,:], color=f'{colours[i]}', 
                label=f'{data_50.experiment[i].values} 50th quantile')
        ax.fill_between(data_50.year, data_90[i,:], data_10[i,:], alpha=0.1, color=f'{colours[i]}', 
                label=f'{data_50.experiment[i].values} 10th and 90th quantile range')

    ax.set_xlim(1985,2050)
    ax.set_title(f'CMIP6 annual {variable} for Brazil')
    ax.set_ylabel(f'{variable}')
    ax.set_xlabel('year')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    ax.grid(linestyle='--')

    fig.savefig(f'{DATADIR}/plots/{variable}.png')