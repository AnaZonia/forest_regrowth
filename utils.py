

# Deleting assets in bulk from projects. Taken from https://gis.stackexchange.com/questions/467363/batch-deleting-of-earth-engine-assets

asset_list = ee.data.listAssets("projects/ee-ana-zonia/assets")["assets"]

def conditional_asset_rm(asset, starts_with):
    """Deletes asset if starts with starts_with """
    id = asset["id"]              # users/username/file  or projects/project-name/assets/file
    findex = 5 if id.startswith("users") else 3
    name = asset["name"]          # projects/earthengine-legacy/assets/users/username/file or projects/project-name/assets/file
    f = name.split("/")[findex]  # file
    if (f.startswith(starts_with)):
        ee.data.deleteAsset(id)
        return f"Deleted asset {id}"

    [conditional_asset_rm(asset, "unif") for asset in asset_list]


def delete_assets_starting_with(starts_with):
    """
    Deletes assets from the given asset list that start with the specified string.

    Args:
        asset_list (list): List of asset dictionaries.
        starts_with (str): The prefix string to match.

    Returns:
        list: List of success messages for deleted assets.
    """
    asset_list = ee.data.listAssets("projects/ee-ana-zonia/assets")["assets"]
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
