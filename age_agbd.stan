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
  // real<lower=0> age_par; // Growth term
}

model {
  // Priors
  // age_par ~ normal(15, 1000);//weibull(1.851616, 111.5);
  // A ~ normal(120, 1000);   // Prior for Asymptote
  // theta ~ normal(0, 1000); // Prior for Shape term
  // sigma ~ exponential(0.01);
  // B0 ~ normal(40, 1000);  // Prior for Intercept
// percentage of the true value is how sigma should be thought of.
  // vector[n] m;
  // m = B0 + A * (1 - exp(-(age_par*age)))^theta;
  agbd ~ normal(B0 + A * (1 - exp(-(0.1*age)))^theta, sigma);
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
