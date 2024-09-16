# scripts/random_forest_fit.py
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold

def fit_random_forest(X, y, n_estimators=200):
    model = RandomForestRegressor(n_estimators=n_estimators, random_state=42)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(model, X, y, cv=kf, scoring='r2')
    
    model.fit(X, y)
    feature_importances = model.feature_importances_
    
    return model, cv_scores, feature_importances
