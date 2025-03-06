import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Optional, List

from xgboost import XGBRegressor
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, GridSearchCV, cross_validate
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder
from sklearn.inspection import permutation_importance
from sklearn.compose import ColumnTransformer

class XGBoostAnalyzer:
    """
    A class for XGBoost regression analysis, including data loading, preprocessing,
    model training, and evaluation.

    Args:
        filepath (str): Path to the data file.
        param_grid (Optional[Dict[str, List]]): Hyperparameter grid for XGBoost.
        feature_columns (Optional[List[str]]): List of features to use.
        biome (int): Biome to analyze (1 for Amazon, 4 for Atlantic Forest).
        final_sample_size (int): Sample size after stratified sampling.
        columns_to_drop (Optional[List[str]]): Columns to drop from the data.
        categorical_features (List[str]): List of categorical feature names.

    Example:
        >>> analyzer = XGBoostAnalyzer(
        ...     filepath="data.csv",
        ...     final_sample_size = 1000,
        ...     feature_columns=["feature1", "feature2"]
        ... )
        >>> analyzer.load_and_preprocess()
        >>> analyzer.train_model()
        >>> analyzer.report_results()
    """

    def __init__(
        self,
        filepath: str,
        param_grid: Optional[Dict] = None,
        feature_columns: Optional[List[str]] = None,
        biome: int = 1,
        final_sample_size: int = 10000,
        cols_to_drop: Optional[List[str]] = None,
        categorical_features: Optional[List[str]] = None
    ):

        self.filepath = filepath
        # xgb__ is necessary since XGBoost is part of a pipeline
        self.param_grid = param_grid or {
            'xgb__learning_rate': [0.01, 0.1],
            'xgb__n_estimators': [100, 500],
            'xgb__max_depth': [3, 7],
            'xgb__min_child_weight': [1, 5]
        }
        self.feature_columns = feature_columns
        self.biome = biome,
        self.final_sample_size = final_sample_size
        self.cols_to_drop = cols_to_drop
        self.categorical_features = categorical_features or ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']
        self.X = None
        self.y = None
        self.results = pd.DataFrame()
        self.feature_importance = None
        self.best_model = None

    def load_and_preprocess(self) -> None:
        """
        Load and clean the data.
        - Remove columns with few non-zero values
            - Ensures extremely rare land use types are not considered
        - Convert categorical features to category type
            - Removes rare categories (fewer than 50 samples)
        - Samples rows to desired sample size if needed.
        """
        
        df = pd.read_csv(self.filepath)
        # Remove columns with yearly climate data (such as "aet_2015" for example)
        df = df.loc[:, ~df.columns.str.contains(r"\d{4}")]
        if self.cols_to_drop:
            df = df.drop(columns = self.cols_to_drop)
        # Keep only rows with the specified biome and remove the 'biome' column
        df = df[df['biome'] == self.biome].drop(columns = 'biome')

        # Sampling
        if self.final_sample_size and len(df) > self.final_sample_size:
            df = df.sample(n = self.final_sample_size, random_state = 42)

        # Identify columns with fewer than 100 non-zero values
        rare_columns = [col for col in df.columns if (df[col] != 0).sum() < 100]
        # Remove columns with few non-zero values
        df = df.drop(columns = rare_columns)
        # Remove rows that had non-zero values in the dropped columns
        df = df[~(df[rare_columns] != 0).any(axis = 1)]

        # Convert categorical features and remove rare categories
        for col in self.categorical_features:
            if col in df.columns:
                df[col] = df[col].astype('category')
                # Identify rare categories (fewer than 50 occurrences)
                counts = df[col].value_counts()
                rare_categories = counts[counts < 50].index
                # Remove rows with rare categories
                df = df[~df[col].isin(rare_categories)]

        # If no parameters are specified, use all columns except 'biomass'
        if self.feature_columns is None:
            self.feature_columns = [col for col in df.columns if col != "biomass"]

        self.X = df[self.feature_columns]
        self.y = df['biomass'].values


    def train_model(self) -> None:
        """
        Train and evaluate the XGBoost model using nested cross-validation.

        This method performs the following steps:
        1. Set up the preprocessing pipeline with MinMaxScaler and OneHotEncoder
        2. Perform nested cross-validation with GridSearchCV
        3. Calculate performance metrics
        4. Compute feature importance using permutation importance

        The best model and results are stored in the class attributes.
        """
        numeric_features = [col for col in self.feature_columns if col not in self.categorical_features]
        categorical_features = [col for col in self.feature_columns if col in self.categorical_features]

        preprocessor = ColumnTransformer(
            transformers=[
                ('num', MinMaxScaler(), numeric_features),
                ('cat', OneHotEncoder(handle_unknown='ignore'), categorical_features)
            ])

        pipeline = Pipeline([
            ('preprocessor', preprocessor),
            ('xgb', XGBRegressor(random_state = 42, use_label_encoder = False, eval_metric='rmse'))
        ])

        inner_kf = KFold(n_splits = 5, shuffle = True, random_state = 42)
        outer_kf = KFold(n_splits = 5, shuffle = True, random_state = 42)

        # Hyperparameter tuning and cross-validation
        grid_search = GridSearchCV(pipeline, self.param_grid, cv = inner_kf,
                                scoring = 'neg_mean_squared_error', refit = 'mse')
        
        cv_results = cross_validate(grid_search, self.X, self.y, cv = outer_kf,
                            scoring = ['neg_mean_squared_error', 'r2'], return_estimator = True)

        # Permutation Importance
        best_model_index = np.argmax(cv_results['test_r2'])
        best_model = cv_results['estimator'][best_model_index]
        self.feature_importance = permutation_importance(best_model, self.X, self.y, n_repeats = 5, random_state = 42)

        # Store results
        self.results = pd.DataFrame([{
            'datasource': self.filepath,
            'r2_mean': np.mean(cv_results['test_r2']),
            'r2_sd': np.std(cv_results['test_r2'])
        }])


    def report_results(self) -> None:
        """
        Visualize and print the results of the analysis.

        This method performs the following steps:
        1. Print regression results (R2 and RMSE)
        2. Print feature importance
        3. Plot feature importance
        """
        # Print numerical results
        print("\nRegression Results:")
        print(self.results.to_string(index = False))
        
        # Show feature importance
        fi_df = pd.DataFrame({
            'feature': self.X.columns,
            'importance': self.feature_importance.importances_mean,
            'std': self.feature_importance.importances_std
        }).sort_values('importance', ascending = False)

        print("\nFeature Importance:")
        print(fi_df.to_string(index = False))

        # Plot feature importance
        plt.figure(figsize = (10, 6))
        plt.barh(fi_df['feature'], fi_df['importance'], 
                xerr = fi_df['std'], color = 'skyblue', edgecolor = 'black')
        plt.xlabel('Importance')
        plt.title('XGBoost Feature Importance')
        plt.gca().invert_yaxis()
        plt.grid(axis = 'x', linestyle = '--', alpha = 0.7)
        plt.tight_layout()
        plt.show()

    def save_results(self, output_path: str) -> None:
        """Save analysis results to CSV."""
        self.results.to_csv(output_path, index = False)


if __name__ == "__main__":
    # Example usage
    analyzer = XGBoostAnalyzer(
        filepath = "~/Documents/data/mapbiomas_heinrich_lulc_regions_renamed.csv",
        final_sample_size = 8000
    )

    analyzer.load_and_preprocess()
    analyzer.train_model()
    analyzer.report_results()
    # analyzer.save_results("./0_results/xgboost_results.csv")