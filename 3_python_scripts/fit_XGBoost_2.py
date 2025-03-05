import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Optional, List, Tuple

from xgboost import XGBRegressor
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, cross_validate, GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance

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

    Example:
        >>> analyzer = XGBoostAnalyzer(
        ...     filepath="data.csv",
        ...     final_sample_size=1000,
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
        cols_to_drop: Optional[List[str]] = None
    ):

        self.filepath = filepath
        self.param_grid = param_grid or {
            'learning_rate': [0.01, 0.1],
            'n_estimators': [100, 500],
            'max_depth': [3, 7],
            'min_child_weight': [1, 5]
        }
        self.feature_columns = feature_columns
        self.biome = biome,
        self.final_sample_size = final_sample_size
        self.cols_to_drop = cols_to_drop
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
            df = df.sample(n=self.final_sample_size, random_state=42)

        # Identify columns with fewer than 100 non-zero values
        rare_columns = [col for col in df.columns if (df[col] != 0).sum() < 100]
        # Remove columns with few non-zero values
        df = df.drop(columns = rare_columns)
        # Remove rows that had non-zero values in the dropped columns
        df = df[~(df[rare_columns] != 0).any(axis = 1)]

        # Convert categorical features and remove rare categories
        for col in ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']:
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
        Train and evaluate the XGBoost model.
        Permutation Importance is used to assess feature importance - measures the increase in the model's prediction error after randomly shuffling the values of a feature.
        """
        pipeline = Pipeline([
            ('scaler', MinMaxScaler()),
            ('xgb', XGBRegressor(random_state = 42))
        ])
        # Define the inner cross-validation strategy (for hyperparameter tuning)
        inner_cv = KFold(n_splits=4, shuffle=True, random_state=42)
        # Define the outer cross-validation strategy (for performance evaluation)
        outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        splits = list(outer_cv.split(self.X))

        # Hyperparameter tuning
        grid_search = GridSearchCV(pipeline, self.param_grid, cv = inner_cv,
                                scoring = 'neg_mean_squared_error', refit = True)

        # Cross-validation
        cv_results = cross_validate(grid_search, self.X, self.y, cv = outer_cv,
                                scoring = ['neg_mean_squared_error', 'r2'],
                                return_estimator = True)
        # mse_scores = -cv_results['test_neg_mean_squared_error']
        r2_scores = cv_results['test_r2']

        # # Permutation importance
        # best_index = np.argmin(-cv_results['test_neg_mean_squared_error'])
        # _, test_indices = splits[best_index]
        # X_test = self.X.iloc[test_indices]
        # X_test = self.X.iloc[kf.split(self.X)[best_index][1]] # Access test indices correctly
        # y_test = self.y[test_indices]
        # best_model = cv_results['estimator'][best_index]
        # self.feature_importance = permutation_importance(best_model, X_test, y_test, n_repeats = 5, random_state = 42)

        # Store results
        self.results = pd.DataFrame([{
            'datasource': self.filepath,
            'r2_mean': np.mean(r2_scores),
            'r2_sd': np.std(r2_scores),
            'r2_final': r2_score(self.y, self._collect_predictions(cv_results, splits))
        }])



    def _collect_predictions(self, cv_results: Dict, splits: List[Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
        """Collect predictions from all cross-validation folds."""
        # # Aggregate predictions from all cross-validation folds.
        # predictions = np.empty(len(self.X))
        # predictions[:] = np.nan
        # for i, estimator in enumerate(cv_results['estimator']):
        #     _, test_indices = list(kf.split(self.X))[i]
        #     fold_predictions = estimator.predict(self.X.iloc[test_indices])
        #     predictions[test_indices] = fold_predictions

        predictions = np.empty(len(self.X))
        predictions[:] = np.nan
        for i, (train_index, test_index) in enumerate(splits):
            fold_predictions = cv_results['estimator'][i].predict(self.X.iloc[test_index])
            predictions[test_index] = fold_predictions
        return predictions



    def report_results(self) -> None:
        """
        Visualize results dataframe (R2 scores)
        Visualize and plot feature importance
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