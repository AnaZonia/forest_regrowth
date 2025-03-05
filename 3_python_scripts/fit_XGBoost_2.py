import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Optional, List

from xgboost import XGBRegressor
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, cross_validate, GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.inspection import permutation_importance
from sklearn.model_selection import StratifiedShuffleSplit


class XGBoostAnalyzer:
    """
    A specialized class for XGBoost regression analysis including data loading,
    preprocessing, model training, and evaluation.
    """

    def __init__(
        self,
        filepath: str,
        param_grid: Optional[Dict] = None,
        pars: Optional[List[str]] = None,
        final_sample_size: int = 10000
    ):
        """
        Initialize the XGBoost analyzer.

        Args:
            filepath (str): Path to the data file
            param_grid (Optional[Dict]): Hyperparameter grid for XGBoost
            pars (Optional[List[str]]): List of features to use
            final_sample_size (int): Sample size after stratified sampling
        """
        self.filepath = filepath
        self.param_grid = param_grid or {
            'learning_rate': [0.01, 0.1],
            'n_estimators': [100, 500],
            'max_depth': [3, 7],
            'min_child_weight': [1, 5]
        }
        self.pars = pars
        self.final_sample_size = final_sample_size
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
        df = df.loc[:, ~df.columns.str.contains(r"\d{4}")] # Removes columns with yearly climate data (such as "aet_2015" for example)

        # Identify columns with fewer than 100 non-zero values
        rare_columns = [col for col in df.columns if (df[col] != 0).sum() < 100]
        # Remove columns with few non-zero values
        df = df.drop(columns = rare_columns)
        # Remove rows that had non-zero values in the dropped columns
        df = df[~(df[rare_columns] != 0).any(axis = 1)]

        # Convert categorical features
        for col in ['topography', 'ecoreg', 'indig', 'protec', 'last_LU']:
            if col in df.columns:
                df[col] = df[col].astype('category')
                # Identify rare categories (fewer than 50 occurrences)
                counts = df[col].value_counts()
                rare_categories = counts[counts < 50].index
                # Remove rows with rare categories
                df = df[~df[col].isin(rare_categories)]

        # Handle feature selection
        if self.pars is None:
            self.pars = [col for col in df.columns if col != "biomass"]

        # Stratified sampling
        if self.final_sample_size and len(df) > self.final_sample_size:
            splitter = StratifiedShuffleSplit(n_splits = 1, 
                test_size = self.final_sample_size, random_state = 42)
            _, sample_index = next(splitter.split(df, df['biome']))
            df = df.iloc[sample_index]

        self.X = df[self.pars]
        self.y = df['biomass'].values


    def train_model(self) -> None:
        """
        Train and evaluate the XGBoost model.
        
        """
        
        pipeline = Pipeline([
            ('scaler', MinMaxScaler()),
            ('xgb', XGBRegressor(random_state = 42))
        ])
        kf = KFold(n_splits = 5, shuffle = True, random_state = 42)
        splits = list(kf.split(self.X))

        # Hyperparameter tuning
        grid_search = GridSearchCV(pipeline, self.param_grid, cv = kf,
                                scoring = 'neg_mean_squared_error', refit = True)
        grid_search.fit(self.X, self.y)
        self.best_model = grid_search.best_estimator_

        # Cross-validation
        cv_results = cross_validate(self.best_model, self.X, self.y, cv = kf,
                                scoring = ['neg_mean_squared_error'],
                                return_estimator = True)
        mse_scores = -cv_results['test_neg_mean_squared_error']
        r2_scores = 1 - (mse_scores / np.var(self.y))

        # Permutation importance
        best_idx = np.argmin(mse_scores)
        X_test = self.X.iloc[splits[best_idx][1]]
        y_test = self.y[splits[best_idx][1]]
        self.feature_importance = permutation_importance(
            self.best_model, X_test, y_test, n_repeats = 10, random_state = 42)

        # Store results
        self.results = pd.DataFrame([{
            'datasource': self.filepath,
            'r2_mean': np.mean(r2_scores),
            'r2_sd': np.std(r2_scores),
            'r2_final': r2_score(self.y, self._collect_predictions(cv_results, splits))
        }])

    def _collect_predictions(self, cv_results: Dict, splits: List) -> np.ndarray:
        """Aggregate predictions from all cross-validation folds."""
        predictions = np.full_like(self.y, np.nan)
        for i, (_, test_idx) in enumerate(splits):
            predictions[test_idx] = cv_results['estimator'][i].predict(self.X.iloc[test_idx])
        return predictions

    def report_results(self) -> None:
        """Print and visualize analysis results."""
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
        self._plot_feature_importance(fi_df)

    def _plot_feature_importance(self, fi_df: pd.DataFrame) -> None:
        """Visualize feature importance."""
        plt.figure(figsize = (10, 6))
        plt.barh(fi_df['feature'], f i_df['importance'], 
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