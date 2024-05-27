r"""
Date
    05/05/2024

Purpose
    This script runs. Hopefully.

Author
    Ana Catarina Avila
    McGill University
"""

# data = import_data("./data/unified_data_15_years.csv")

# old_data = pd.read_csv("./data/amazon_df_sample_10million.csv")

# # Define the years and climatic variables
# years = np.arange(1985, 2020)

# # define the climatic parameters - the ones that change yearly
# climatic = ["prec", "si"]
# # define the non-climatic parameters - the ones that are fixed throughout regrowth and
# # that are used for fitting the model (excludes age and agbd)
# non_climatic = [
#     col for col in data.columns if not any(var in col for var in climatic + ["agbd"])
# ]



# # write the growth curve with yearly climatic data and permanent non-climatic data
# def growth_curve(pars, data, years):
#     r"""This function defines the growth function and parameter dictionary"""

#     # # Calculate k
#     # k = np.zeros(len(data))
#     # for year in years:
#     #     for clim_var in climatic:
#     #         k += pars[clim_var] * data[f"{clim_var}_{year}"]
#     #     for unique_var in non_climatic:
#     #         k += pars[unique_var] * data[unique_var]

#     return (
#         pars["B0"]
#         + pars["A"]
#         * (1 - np.exp(-pars["last_fire"] * data["last_fire"])) ** pars["theta"]
#     )


# def negative_log_likelihood(pars, data):
#     # Assuming a normal distribution for the data
#     likelihood = norm.pdf(
#         data["agbd"] - growth_curve(pars, data, years), loc=0, scale=pars["sd"]
#     )
#     return -np.sum(np.log(likelihood))



def main():
    path = "./data/unified_data_15_years.csv"
    import_data(path)


#     years = np.arange(1985, 2020)
#     # Initialize parameters
#     pars = {'B0': 40, 'A': 80, 'age': 0.05, 'theta': 1.5}

#     run_model()

    
    # x0 = list(pars.values())


    # # Now you can call minimize with the modified growth_curve function
    # result = minimize(negative_log_likelihood, x0, args=(data,), method="Nelder-Mead")


if __name__ == "__main__":
    main()
