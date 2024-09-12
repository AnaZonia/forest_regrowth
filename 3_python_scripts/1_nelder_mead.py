# scripts/nelder_mead_fit.py
from scipy.optimize import minimize
import numpy as np
import pandas as pd

import pandas as pd
import numpy as np
from typing import List



# write the growth curve with yearly climatic data and permanent non-climatic data
def growth_curve(pars, X, y):

    # Extract parameters excluding 'B0' and 'theta' as a DataFrame
    data_pars = pars.drop(columns=['B0', 'theta'], errors='ignore')
    print(data_pars)

    # Calculate k using vectorized operations
    k = np.dot(data[data_pars.index].values, data_pars.values) + 1

    # Calculate and return the growth curve
    B0 = pars['B0']
    theta = pars['theta']
    nearest_mature = data['nearest_mature'].values
    return B0 + (nearest_mature - B0) * (1 - np.exp(-k))**theta




def run_optim(train_data, pars, conditions, test_data=None):
    if 'age' in pars.index:
        conditions.append(lambda pars: pars['age'] < 0)
    
    def objective(params):
        return likelihood(pd.Series(params, index=pars.index), train_data, conditions)
    
    result = minimize(objective, pars.values, method='Nelder-Mead')
    
    return result


def find_combination_pars(pars: List[str], data: pd.DataFrame) -> pd.Series:
    # Initialize parameter DataFrame with basic parameters and theta
    all_pars_iter = pd.DataFrame(0, index=[0], columns=pars)
    all_pars_iter["theta"] = 1
    all_pars_iter["B0"] = data["agbd"].mean()

    # Initialize the best model with basic parameters
    remaining = list(range(len(pars)))
    taken = len(remaining) + 1  # out of range of values such that remaining[-taken:] = [] for the first iteration

    # best model dictionary
    best = {"AIC": 0}

    basic_pars = ["theta", "B0"]
    # beginning with the essential parameters
    best["par"] = all_pars_iter[basic_pars].iloc[0]
    val = len(basic_pars)

    base_row = pd.DataFrame({"likelihood": 0}, index=[0])
    base_row = pd.concat([base_row, all_pars_iter], axis=1)
    base_row[:] = np.nan

    # Iteratively add parameters and evaluate the model. Keep only AIC improvements.
    for i in range(len(pars)):
        iter_df = pd.DataFrame()
        for j in remaining[:-taken]:
            print(j)
            inipar = pd.concat([best["par"], all_pars_iter[pars[j]]])

            model = run_optim(data, inipar, conditions)

            iter_row = base_row.copy()
            iter_row.loc[0, model["par"].index] = model["par"]
            iter_row.loc[0, "likelihood"] = model["value"]
            iter_df = pd.concat([iter_df, iter_row], ignore_index=True)


        best_model = iter_df["likelihood"].idxmin()
        best_model_AIC = 2 * iter_df.loc[best_model, "likelihood"] + 2 * (i + val + 1)

        if best["AIC"] == 0 or best_model_AIC < best["AIC"]:
            best["AIC"] = best_model_AIC
            best["par"] = iter_df.loc[best_model, all_pars_iter.columns].dropna()
            taken = [i for i, par in enumerate(pars) if par in best["par"].index]
        else:
            print("No improvement. Exiting loop.")
            break

    return best["par"]
