data {
  int<lower=0> n; // Number of observations
  real ages[n];
  real agbds[n];
}

parameters {
//  vector[M] beta; // Coefficients for predictor variables
  real<lower=0, upper = 100> B0; // Intercept
  real<lower=0, upper = 400> A;  // Asymptote
  real<lower=0> theta; // Shape term
  real<lower=0> age; // Growth term
  real<lower=0> sigma; // process error
}

model {
  // Priors
  age ~ normal(0, 1); // Assuming standard normal priors for beta coefficients
  B0 ~ normal(0, 5);  // Prior for Intercept
  A ~ normal(0, 10);   // Prior for Asymptote
  theta ~ normal(0, 10); // Prior for Shape term
  sigma ~ cauchy(0,2.5);

  // Likelihood
  for (pixel in 1:n){
    agbds[pixel] ~ normal(B0 + A * (1 - exp(-(age*ages[pixel])))^theta, sigma);
  }
}
