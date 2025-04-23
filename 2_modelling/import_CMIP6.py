
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
                    "month": [f"{m:02d}" for m in range(1, 13)]  # All months
                    # "area": [5.3, -74, -35, -34]
                }).download(f'{DATADIR}/{variable}/{variable}_{experiment}_{start}-{end-1}_{model}.zip')
                
                print(f"✅ Downloaded {variable} for {model} ({experiment}, {start}-{end-1})")

            except Exception as e:
                print(f"⚠️ Failed to download {variable} for {model} ({experiment}, {start}-{end-1}): {e}")

def unzip_all(variable):
    """Unzip all downloaded zip files for a variable."""
    zip_paths = glob(f'{DATADIR}/{variable}/*.zip')
    for zp in zip_paths:
        with zipfile.ZipFile(zp, 'r') as zip_ref:
            zip_ref.extractall(f'{DATADIR}/{variable}/')
    print(f"✅ Extracted all zip files for {variable}")

# Function to merge NetCDF files
def merge_ncdf(files, output_path):
    """Merge multiple NetCDF files along time dimension and save."""
    datasets = [xr.open_dataset(f) for f in files]
    combined = xr.concat(datasets, dim="time").sortby('time')  # Merge along time dimension
    combined.to_netcdf(output_path)
    print(f"Saved merged dataset: {output_path}")

def merge_all(variable):
    """Find NetCDF files for a variable, group by model/experiment, and merge."""
    data_files = sorted(glob(f"{DATADIR}/{variable}/*.nc"))
    grouped_files = {}
    for file in data_files:
        base_name = "_".join(os.path.basename(file).split('_')[:-1])  # Remove year part
        grouped_files.setdefault(base_name, []).append(file)
    for base_name, files in grouped_files.items():
        output_file = f"{DATADIR}/{variable}/{base_name}_merged.nc"
        merge_ncdf(files, output_file)



def geog_agg(nc_file, variable_short):
    """Aggregate data spatially (weighted mean) and temporally (annual mean), convert units."""
    ds = xr.open_dataset(nc_file, engine='netcdf4')
    exp = ds.attrs.get('experiment_id', 'unknown')
    mod = ds.attrs.get('source_id', 'unknown')
    da = ds[variable_short]
    weights = np.cos(np.deg2rad(da.lat))
    da_weighted = da.weighted(weights)
    da_agg = da_weighted.mean(['lat', 'lon'])
    da_yr = da_agg.groupby('time.year').mean()
    # Convert Kelvin to Celsius if temperature
    if variable_short in ['ta', 'tas', 'air_temperature', 'daily_maximum_near_surface_air_temperature']:
        da_yr = da_yr - 273.15
    da_yr = da_yr.assign_coords(model=mod)
    da_yr = da_yr.expand_dims('model')
    da_yr = da_yr.assign_coords(experiment=exp)
    da_yr = da_yr.expand_dims('experiment')
    out_path = f'{DATADIR}/{variable_short}/agg_{exp}_{mod}_{str(da_yr.year[0].values)}.nc'
    da_yr.to_netcdf(path=out_path)
    print(f"Saved aggregated data: {out_path}")






def plot_ensemble(variable, experiments, variable_short):
    """Plot ensemble median and quantiles of annual global mean variable."""

    data_ds = xr.open_mfdataset(f'{DATADIR}/{variable}/agg*.nc')
    data_ds.load()
    data = data_ds[f'{variable_short}']

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

    ax.set_xlim(1985, 2050)
    ax.set_title(f'CMIP6 annual global {variable_short}')
    ax.set_ylabel(f'{variable_short} (°C)')
    ax.set_xlabel('Year')
    ax.legend()
    ax.grid(linestyle='--')
    plot_path = f'{DATADIR}/plots/{variable_short}.png'
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    fig.savefig(plot_path)
    print(f"Plot saved to {plot_path}")



def save_raster_summary(variable_short, variable):
    """Save a raster dataset summarizing the ensemble median across models and experiments."""
    agg_files = glob(f'{DATADIR}/{variable}/*merged.nc')
    if not agg_files:
        print("No merged files to create raster summary.")
        return
    # Load all merged datasets
    ds_all = xr.open_mfdataset(agg_files, combine='by_coords')
    da = ds_all[variable_short]
    # Convert Kelvin to Celsius if temperature
    if variable_short in ['ta', 'tas', 'air_temperature', 'daily_maximum_near_surface_air_temperature']:
        da = da - 273.15
    # Compute median across models and experiments for each time and grid cell
    da_median = da.median(dim=['model', 'experiment'])
    # Save raster summary as NetCDF
    raster_path = f'{DATADIR}/{variable_short}/ensemble_median_raster.nc'
    da_median.to_netcdf(raster_path)
    print(f"Saved ensemble median raster data to {raster_path}")

def main(download=False, plot=False, save_raster=False):
    """Main workflow to download, process, plot, and save raster data."""
    for variable_short, variable in variables.items():
        os.makedirs(f'{DATADIR}/{variable}', exist_ok=True)
        if download:
            # Download historical data
            download_data("historical", 1985, 2015, models, variable)
            # Download future scenarios
            for exp in experiments[1:]:
                download_data(exp, 2015, 2100, models, variable)
            unzip_all(variable)
            merge_all(variable)
        aggregate_all(variable_short)
        if plot:
            plot_ensemble(variable_short, experiments)
        if save_raster:
            save_raster_summary(variable_short)

# Load and prepare data
def save_nc(variables, experiments, models):
    for variable_short, variable in variables.items():
        # DOWNLOAD DATA FOR HISTORICAL PERIOD
        download_data("historical", 1985, 2015, models, variable)

        # DOWNLOAD DATA FOR FUTURE SCENARIOS
        for i in experiments[1:]:
            download_data(i, 2015, 2100, models, variable)

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


        # Merge files and save
        for base_name, files in grouped_files.items():
            output_file = f"{DATADIR}/{variable}/{base_name}_merged.nc"
            merge_ncdf(files, output_file)

        cmip6_nc = list()
        cmip6_nc_rel = glob(f'{DATADIR}/{variable}/*merged.nc')
        for i in cmip6_nc_rel:
            cmip6_nc.append(os.path.basename(i))


if __name__ == "__main__":
    # # Download and process data
    # save_nc(variables, experiments, models)

    # # Plot data
    # for variable_short, variable in variables.items():
    #     plot_data(variable, experiments, variable_short)

    # Save raster summary
    for variable_short, variable in variables.items():
        save_raster_summary(variable_short, variable)