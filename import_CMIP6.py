
# Tutorial at https://ecmwf-projects.github.io/copernicus-training-c3s/projections-cmip6.html

# CMIP6 climate projections https://cds.climate.copernicus.eu/datasets/projections-cmip6?tab=download
# Agroclimatic indicators derived https://cds.climate.copernicus.eu/datasets/sis-agroclimatic-indicators?tab=download

# note that the request may be queued, so it may take a while to download
# if in doubt, check in the Copernicus Climate Data Store (CDS) for the status of your request

# General libs for file paths\, data extraction, etc
import urllib3 
import cdsapi
from pathlib import Path

client = cdsapi.Client()
urllib3.disable_warnings() # Disable warnings for data download via API

DATADIR = "./0_data/CMIP6"

experiments = ['historical', 'ssp1_2_6', 'ssp2_4_5', 'ssp5_8_5']

models = ['hadgem3_gc31_ll', 'inm_cm5_0', 'inm_cm4_8', 'ipsl_cm6a_lr', 
          'miroc_es2l', 'mpi_esm1_2_lr', 'ukesm1_0_ll']

variables = {
    "ta": "air_temperature",
    "pr": "precipitation",
    "rsds": "surface_downwelling_shortwave_radiation",
    "tas": "near_surface_air_temperature",
    "huss": "near_surface_specific_humidity",
    "mrsos": "moisture_in_upper_portion_of_soil_column"
}

def download_data(experiment, start_year, end_year, models, variable = "air_temperature"):
    """Download CMIP6 data from CDS for given experiment, years, models and variable."""

    # Create a directory for the variable if it doesn't exist
    var_dir = Path(DATADIR) / variable
    var_dir.mkdir(parents=True, exist_ok=True)

    for model in models:
        for start in range(start_year, end_year, 10):
            end = min(start + 10, end_year)  # Ensure no year exceeds end_year
            year_range = [str(y) for y in range(start, end)]
            
            filename = f"{variable}_{experiment}_{start}-{end}_{model}.zip"
            filepath = var_dir / filename

            if filepath.is_file():
                print(f"ℹ️ File already exists: {filepath}, skipping download.")
                continue

            request = {
                    "temporal_resolution": "monthly",
                    "experiment": experiment,
                    "variable": variable,
                    "model": model,
                    "year": year_range,
                    "month": [f"{m:02d}" for m in range(1, 13)],  # All months
                    "area": [5.3, -74, -35, -34]
                }

            # Add level only for specific variables
            if variable in ["air_temperature"]:
                request["level"] = ["1000"]

            try:
                client.retrieve("projections-cmip6", request).download(str(filepath))
                print(f"✅ Downloaded {variable} for {model} ({experiment}, {start}-{end})")

            except Exception as e:
                print(f"⚠️ Failed to download {variable} for {model} ({experiment}, {start}-{end}): {e}")


def main():
    """Main function to organize downloads by variable and experiment"""

    for variable_short, variable in variables.items():
        print(f"\n{'='*40}")
        print(f"PROCESSING VARIABLE: {variable_short.upper()} ({variable})")
        print(f"{'='*40}")

        for exp in experiments:
            # Set temporal ranges
            if exp == 'historical':
                start, end = 1950, 2014  # Historical period
            else:
                start, end = 2015, 2075  # Future scenarios

            print(f"\n📥 Downloading {variable} ({exp.upper()} {start}-{end})")
            download_data(
                experiment = exp,
                start_year = start,
                end_year = end,
                models = models,
                variable = variable
            )


if __name__ == "__main__":
    main()