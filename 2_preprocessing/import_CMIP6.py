import urllib3
import cdsapi
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

urllib3.disable_warnings()

client = cdsapi.Client()
DATADIR = "./0_data/CMIP6"

experiments = ['historical']  # Add 'ssp1_2_6', etc. as needed

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


def make_request(experiment, start, end, model, variable):
    """Submit a single data request to CDS."""
    var_dir = Path(DATADIR) / variable
    var_dir.mkdir(parents=True, exist_ok=True)

    year_range = [str(y) for y in range(start, end)]
    filename = f"{variable}_{experiment}_{start}-{end}_{model}.zip"
    filepath = var_dir / filename

    if filepath.is_file():
        print(f"ℹ️ File already exists: {filepath}, skipping download.")
        return

    request = {
        "temporal_resolution": "monthly",
        "experiment": experiment,
        "variable": variable,
        "model": model,
        "year": year_range,
        "month": [f"{m:02d}" for m in range(1, 13)],
        "area": [5.3, -74, -35, -34]
    }

    if variable == "air_temperature":
        request["level"] = ["1000"]

    try:
        client.retrieve("projections-cmip6", request).download(str(filepath))
        print(f"✅ Downloaded {variable} for {model} ({experiment}, {start}-{end})")
    except Exception as e:
        print(f"⚠️ Failed to download {variable} for {model} ({experiment}, {start}-{end}): {e}")


def main():
    tasks = []

    with ThreadPoolExecutor(max_workers=5) as executor:  # Adjust workers as needed
        for variable_short, variable in variables.items():
            print(f"\n{'='*40}")
            print(f"PROCESSING VARIABLE: {variable_short.upper()} ({variable})")
            print(f"{'='*40}")

            for exp in experiments:
                start_year, end_year = (2010, 2015) if exp == "historical" else (2015, 2076)

                for model in models:
                    for start in range(start_year, end_year, 10):
                        end = min(start + 10, end_year)
                        tasks.append(executor.submit(
                            make_request, exp, start, end, model, variable
                        ))

        # Optional: monitor task completion
        for future in as_completed(tasks):
            _ = future.result()  # Triggers exceptions if any

if __name__ == "__main__":
    main()
