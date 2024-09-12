# scripts/hyperparameter_tuning.py
from skopt import BayesSearchCV
from skopt.space import Real, Integer
from sklearn.ensemble import RandomForestRegressor
import pandas as pd

def tune_random_forest(X:pd.DataFrame, y):
    param_space = {
        'n_estimators': Integer(100, 1000),
        'max_depth': Integer(3, 20),
        'min_samples_split': Integer(2, 20),
        'min_samples_leaf': Integer(1, 10),
    }
    
    rf = RandomForestRegressor(random_state=42)
    bayes_cv = BayesSearchCV(rf, param_space, n_iter=50, cv=5, 
                             scoring='r2', random_state=42, n_jobs=-1)
    bayes_cv.fit(X, y)
    
    return bayes_cv.best_estimator_, bayes_cv.best_params_, bayes_cv.best_score_


def test_tune_random_forest():
    X = np.random.rand(100, 10)
    y = np.random.rand(100)
    
    model, params, score = tune_random_forest(X, y)
    
    assert isinstance(model, RandomForestRegressor), "Model is not an instance of RandomForestRegressor"
    assert isinstance(params, dict)
    assert isinstance(score, float)
    
    assert len(params) > 0
    assert score > 0.0
    assert score <= 1.0


def run_tune_random_forest():
    # run with real growj up data
    X, y = load_and_preprocess_data()
    model, params, score = tune_random_forest(X, y)
    
    print(f"Best model parameters: {params}, Best R2 score: {score}")
    return model, params, score