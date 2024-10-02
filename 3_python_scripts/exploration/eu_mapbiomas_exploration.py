import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

# Load the CSV file into a DataFrame
df = pd.read_csv("./0_data/mapbiomas_eu.csv")
# Remove the first and last columns
df = df.iloc[:, 1:-1]

print(df[['age_mapbiomas', 'age_eu']].corr())


# Plot median agbd per value of age_eu
median_agbd_age_eu = df.groupby('age_eu')['agbd'].median().reset_index()
plt.figure(figsize=(10, 5))
sns.lineplot(x='age_eu', y='agbd', data=median_agbd_age_eu, marker='o')
plt.title('Median AGBD per Age EU')
plt.xlabel('age_eu')
plt.ylabel('Median AGBD')
plt.show()

# Plot median agbd per value of age_mapbiomas
median_agbd_age_mapbiomas = df.groupby('age_mapbiomas')['agbd'].median().reset_index()
plt.figure(figsize=(10, 5))
sns.lineplot(x='age_mapbiomas', y='agbd', data=median_agbd_age_mapbiomas, marker='o')
plt.title('Median AGBD per Age Mapbiomas')
plt.xlabel('age_mapbiomas')
plt.ylabel('Median AGBD')
plt.show()

# Perform linear regression with age_eu against agbd
X_age_eu = df[['age_eu']]
y = df['agbd']

# Add a constant to the independent variables
X_age_eu = sm.add_constant(X_age_eu)

# Fit the linear regression model
model_age_eu = sm.OLS(y, X_age_eu).fit()

# Print the summary of the linear regression model
print("Linear Regression Model for age_eu against agbd")
print(model_age_eu.summary())

# Perform linear regression with age_mapbiomas against agbd
X_age_mapbiomas = df[['age_mapbiomas']]

# Add a constant to the independent variables
X_age_mapbiomas = sm.add_constant(X_age_mapbiomas)

# Fit the linear regression model
model_age_mapbiomas = sm.OLS(y, X_age_mapbiomas).fit()

# Print the summary of the linear regression model
print("Linear Regression Model for age_mapbiomas against agbd")
print(model_age_mapbiomas.summary())

