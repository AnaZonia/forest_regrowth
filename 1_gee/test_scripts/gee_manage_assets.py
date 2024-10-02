
import ee
import geemap
from gee_0_utils import *
initialize()

# Deleting assets in bulk from projects. Taken from https://gis.stackexchange.com/questions/467363/batch-deleting-of-earth-engine-assets

asset_list = ee.data.listAssets("projects/amazon-forest-regrowth/assets")["assets"]

def conditional_asset_rm(asset_list, starts_with):
    """Deletes assets from a list if they start with starts_with."""
    success_messages = []
    for asset in asset_list:
        id = asset["id"]
        name = asset["name"]
        findex = 5 if id.startswith("users") else 3
        f = name.split("/")[findex]
        if f.startswith(starts_with):
            ee.data.deleteAsset(id)
            success_messages.append(f"Deleted asset {id}")
    return success_messages


def move_assets_to_folder(starts_with, destination_folder):
    """
    Moves assets from the specified asset list to a given folder 
    if their names start with the specified prefix.

    Args:
        starts_with (str): The prefix string to match.
        destination_folder (str): The folder path to move the assets into.

    Returns:
        list: List of success messages for moved assets.
    """
    success_messages = []
    
    for asset in asset_list:
        id = asset["id"]
        name = asset["name"]
        findex = 5 if id.startswith("users") else 3
        f = name.split("/")[findex]
        
        if f.startswith(starts_with):
            # Define the new asset ID within the destination folder
            new_id = f"{destination_folder}/{f}"
            # Move the asset to the new ID
            ee.data.copyAsset(id, new_id)
            ee.data.deleteAsset(id)
            success_messages.append(f"Moved asset {id} to {new_id}")
    
    return success_messages


# move_assets_to_folder("near", "projects/amazon-forest-regrowth/assets/mapbiomas")

# conditional_asset_rm(asset_list, "land")