data {
  int <lower=0>n;// Number of observations
  vector[n] age;// Corrected with semicolon
  vector[n] agbd;// Corrected with semicolon
}

parameters {
  real<lower=0> B0; // Intercept
  real<lower=0> A;  // Asymptote
  real<lower=0> theta; // Shape term
  real<lower=0> sigma; // process error
  real<lower=0> k; // Growth term
}

transformed parameters {
    real<lower=0, upper=1> p;
    p = B0 / A;

}

model {
  // Priors
  k ~ normal(0.001, 1);
  A ~ normal(160, 60);   // Prior for Asymptote
  theta ~ normal(5, 3); // Prior for Shape term
  sigma ~ exponential(0.1);
  // B0 ~ normal(40, 10);  // Prior for Intercept
  p ~ beta(3,5);

  // vector[n] m;
  // m = B0 + A * (1 - exp(-(age_par*age)))^theta;
  agbd ~ normal(A * (p + (1 - exp(-(k*age)))^theta), sigma);

  // agbd ~ log_A + log_sum_exp(log_inv_logit(p), 1/log_inv_logit(k*age))

  // agbd ~ normal(B0 * (1 + (A/B0) * (1 - exp(-(k*age)))^theta), sigma);

// either set parameters with a parameter with a lower bound of 
// or add a value with a distribution ranging from zero to the maximum year of colonization
}

// B0 + exp(log_A) * (1 - exp(-(age_par*age)))^theta

// generated quantities {
//   vector[n] agbd; // Corrected with semicolon
//   vector[n] mean_agbd; // Corrected with semicolon
//   mean_agbd = B0 + exp(log_A) * (1 - exp(-(age_par*age)))^theta;
//   for (pixel in 1:n){
//     agbd[pixel] = normal_rng(mean_agbd[pixel],sigma);
//   }
// }
