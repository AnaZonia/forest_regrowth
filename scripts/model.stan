data {
  int<lower=1> N; // Number of observations
  int<lower=1> M; // Number of predictor variables (climatic and non-climatic)
  matrix[N, M] X; // Predictor matrix
  vector[N] agbd;    // Response variable
}

parameters {
  vector[M] beta; // Coefficients for predictor variables
  real<lower=0, upper = 100> B0; // Intercept
  real<lower=0, upper = 400> A;  // Asymptote
  real<lower=0> theta; // Shape term
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
