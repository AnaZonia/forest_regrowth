
# - Imports the dataframe
# - Removes unnecessary columns that will not be used in analysis
# - Converts categorical data to dummy variables

def import_data(path: str):
    r"""This function:
    - Imports the dataframe
    - Removes unnecessary columns that will not be used in analysis
    - Converts categorical data to dummy variables"""

    # Read the CSV file
    data = pd.read_csv("./data/unified_data_15_years.csv")

    # Drop unnecessary columns
    data = data.drop(columns=["system:index", ".geo", "biome"])

    # Convert 'soil' and 'ecoreg' to categorical data
    categorical = ["soil", "ecoreg"]
    data[categorical] = data[categorical].astype("category")

    # Create dummy variables for 'ecoreg' and 'soil'
    data = pd.get_dummies(data, columns=["ecoreg", "soil"])

    return data