import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from data_utils import load_and_preprocess_data


# Example with increased complexity
class CurveFitModel(nn.Module):
    def __init__(self, input_dim):
        super(CurveFitModel, self).__init__()
        self.hidden1 = nn.Linear(input_dim, 128)
        self.hidden2 = nn.Linear(128, 64)
        self.hidden3 = nn.Linear(64, 32)
        self.output = nn.Linear(32, 1)  # No activation in output layer for regression
        
    def forward(self, x):
        x = torch.relu(self.hidden1(x))
        x = torch.relu(self.hidden2(x))
        x = torch.relu(self.hidden3(x))
        return self.output(x)  # Linear activation for regression


def generate_curve_data(model, X_scaler, y_scaler, feature_index, num_points=100):
    # Create a range of values for the chosen feature
    feature_min, feature_max = X_scaler.data_min_[feature_index], X_scaler.data_max_[feature_index]
    feature_range = np.linspace(feature_min, feature_max, num_points).reshape(-1, 1)
    
    # Create input data for prediction
    X_curve = np.zeros((num_points, X_scaler.n_features_in_))
    X_curve[:, feature_index] = feature_range.flatten()
    
    # Scale the input data
    X_curve_scaled = X_scaler.transform(X_curve)
    
    # Generate predictions
    model.eval()
    with torch.no_grad():
        y_pred_scaled = model(torch.FloatTensor(X_curve_scaled)).numpy()
    
    # Inverse transform the predictions
    y_pred = y_scaler.inverse_transform(y_pred_scaled)
    
    return feature_range, y_pred



def train_curve_fit(X, y, model, epochs=1000, lr=0.01):
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    X_tensor = torch.FloatTensor(X)
    y_tensor = torch.FloatTensor(y).reshape(-1, 1)
    
    for epoch in range(epochs):
        model.train()  # Ensure the model is in training mode
        outputs = model(X_tensor)
        loss = criterion(outputs, y_tensor)
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        if (epoch+1) % 100 == 0:
            print(f'Epoch [{epoch+1}/{epochs}], Loss: {loss.item():.4f}')
    
    return model

def evaluate_curve_fit(model, X, y):
    model.eval()  # Ensure the model is in evaluation mode
    with torch.no_grad():
        X_tensor = torch.FloatTensor(X)
        y_pred = model(X_tensor).detach().cpu().numpy().flatten()  # Ensure it's on CPU
    r2 = r2_score(y, y_pred)
    return r2, y_pred

def curve_fit_main():
    # Define biomes and data sources
    biomes = [1, 4, "both"]
    filepaths = {
        "mapbiomas": "./0_data/non_aggregated.csv"
    }

    results_list = []
    
    for biome in biomes:
        for datasource, filepath in filepaths.items():
            # Load and preprocess data
            X, y, _, _ = load_and_preprocess_data(filepath,
                                biome=biome,
                                final_sample_size=30000,
                                unseen_portion=0.2)
            
            # Split the data
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
            
            # Scale the features
            scaler_X = MinMaxScaler()
            X_train_scaled = scaler_X.fit_transform(X_train)
            X_test_scaled = scaler_X.transform(X_test)
            
            # Scale the target variable
            scaler_y = MinMaxScaler()
            y_train_scaled = scaler_y.fit_transform(y_train.reshape(-1, 1)).flatten()
            y_test_scaled = scaler_y.transform(y_test.reshape(-1, 1)).flatten()
            
            # Create and train the model
            model = CurveFitModel(input_dim=X_train.shape[1])
            trained_model = train_curve_fit(X_train_scaled, y_train_scaled, model, epochs=2000, lr=0.001)

            
            # Evaluate the model on the training set
            train_r2, train_pred_scaled = evaluate_curve_fit(trained_model, X_train_scaled, y_train_scaled)
            # Reverse scaling for predictions
            train_pred = scaler_y.inverse_transform(train_pred_scaled.reshape(-1, 1)).flatten()
            
            # Evaluate the model on the test set
            test_r2, test_pred_scaled = evaluate_curve_fit(trained_model, X_test_scaled, y_test_scaled)
            # Reverse scaling for predictions
            test_pred = scaler_y.inverse_transform(test_pred_scaled.reshape(-1, 1)).flatten()
            
            print(f"\nNeural Network Curve Fit Results for Biome {biome}, DataSource: {datasource}")
            print(f"Train R2: {train_r2:.3f}")
            print(f"Test R2: {test_r2:.3f}")

            feature_index = 0
            feature_name = X.columns[feature_index]
            feature_range, curve_pred = generate_curve_data(trained_model, scaler_X, scaler_y, feature_index)
            # Visualize the fitted curve
            plt.figure(figsize=(12, 6))

            # Plot the actual data points
            plt.scatter(X_train.iloc[:, feature_index], y_train, alpha=0.5, label='Actual Data')  # Fixed here

            # Plot the fitted curve
            plt.plot(feature_range, curve_pred, 'r-', label='Fitted Curve')

            plt.xlabel(feature_name)
            plt.ylabel('Target Variable')
            plt.title(f'Fitted Curve for {feature_name} - Biome {biome} - {datasource}')
            plt.legend()
            plt.show()
           
            # # Visualize the curve fit for test set
            # plt.figure(figsize=(10, 6))
            # plt.scatter(y_test, test_pred, alpha=0.5)
            # plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
            # plt.xlabel('Actual Values')
            # plt.ylabel('Predicted Values')
            # plt.title(f'Neural Network Curve Fit - Biome {biome} - {datasource}')
            # plt.show()
            # Choose a feature to visualize (e.g., the first feature)


if __name__ == "__main__":
    curve_fit_main()


