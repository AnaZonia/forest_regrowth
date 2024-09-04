# main.py

def main():
    X, y = load_and_preprocess_data("data/processed_amaz_10yr.csv")
    
    # Fit Nelder-Mead
    nelder_mead_params = fit_nelder_mead(X, y)
    print(f"Nelder-Mead fitted parameters: {nelder_mead_params}")

    # Fit Random Forest
    rf_model, cv_scores, feature_importances = fit_random_forest(X, y)
    print(f"Random Forest CV scores: {cv_scores}")
    print(f"Random Forest feature importances: {feature_importances}")

    # Hyperparameter tuning
    best_model, best_params, best_score = bayesian_optimization(X, y)
    print(f"Best model parameters: {best_params}, Best R2 score: {best_score}")

    # Evaluate the best model
    mse, r2 = evaluate_model(best_model, X, y)
    print(f"Evaluation - MSE: {mse}, R2: {r2}")

if __name__ == "__main__":
    main()
