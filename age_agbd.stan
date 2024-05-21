data {
  int <lower=0>n;// Number of observations
  real age[n];// Corrected with semicolon
  real agbd[n]; // Corrected with semicolon
}


parameters {
  real<lower=0, upper = 100> B0; // Intercept
  real<lower=100, upper = 400> A;  // Asymptote
  real<lower=0> theta; // Shape term
  real<lower=0> sigma; // process error
  real<lower=0> age_par; // Growth term
}

model {
  // Priors
  age_par ~ normal(0, 1); // Assuming standard normal priors for beta coefficients
  B0 ~ normal(0, 5);  // Prior for Intercept
  A ~ normal(0, 10);   // Prior for Asymptote
  theta ~ normal(0, 5); // Prior for Shape term
  sigma ~ cauchy(0,10);

  // Likelihood
  for (pixel in 1:n){
    agbd[pixel] ~ normal(B0 + A * (1 - exp(-(age_par*age[pixel])))^theta, sigma);
  }
}
