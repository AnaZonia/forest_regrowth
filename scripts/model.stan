data {
  int<lower=1> N; // Number of observations
  int<lower=1> M; // Number of predictor variables (climatic and non-climatic)
  matrix[N, M] X; // Predictor matrix
  vector[N] agbd; // Response variable
  vector[33] years; // Vector of years
  vector[33] climatic; // Vector of climatic variables
  vector[33] non_climatic; // Vector of non-climatic variables
}

parameters {
  vector[M] beta; // Coefficients for predictor variables
  real<lower=0, upper = 100> B0; // Intercept
  real<lower=0, upper = 400> A;  // Asymptote
  real<lower=0> theta; // Shape term
}

functions {
 real growth_curve(vector pars, matrix data, vector years, vector climatic, vector non_climatic) {
    real k = 0;
    for (i in 1:size(years)) {
      for (j in 1:size(climatic)) {
        k = k + pars[climatic[j]] * data[i, climatic[j]];
      }
      for (j in 1:size(non_climatic)) {
        k = k + pars[non_climatic[j]] * data[i, non_climatic[j]];
      }
    }
    return pars[1] + pars[2] * (1 - exp(-k))^pars[3];
 }
}


model {
  // Priors
  beta ~ normal(0, 1); // Assuming standard normal priors for beta coefficients
  B0 ~ normal(0, 50);  // Prior for Intercept
  A ~ normal(0, 100);   // Prior for Asymptote
  theta ~ normal(0, 10); // Prior for Shape term
  
  
  // Likelihood
  agbd ~ normal(B0 + X * beta * (1 - exp(-(X * beta))), A * (1 - exp(-(X * beta)))^theta);
}
