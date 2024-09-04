# scripts/hyperparameter_tuning.py
from skopt import BayesSearchCV
from skopt.space import Real, Integer
from sklearn.ensemble import RandomForestRegressor

def bayesian_optimization(X, y):
    param_space = {
        'n_estimators': Integer(100, 1000),
        'max_depth': Integer(3, 10),
        'min_samples_split': Integer(2, 20),
        'min_samples_leaf': Integer(1, 10),
    }
    
    rf = RandomForestRegressor(random_state=42)
    bayes_cv = BayesSearchCV(rf, param_space, n_iter=50, cv=5, scoring='r2', random_state=42, n_jobs=-1)
    bayes_cv.fit(X, y)
    
    return bayes_cv.best_estimator_, bayes_cv.best_params_, bayes_cv.best_score_
