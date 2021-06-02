//Stan normal model with uncertainty NO RAPI with species interaction
data {
  int<lower=0> N;       // number of cases
  int<lower=1> K;       // number of species
  real y[N];          // outcome (variate)
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  matrix[N,K] x_species; //all species data matrix
  //vector[N] x_Anbo;       // species
  //vector[N] x_Rhma;       // species
  //vector[N] x_Osse;       // species
  //vector[N] x_Raca;       // species

}
parameters {
  vector[K] b_species;
  //real b_Anbo;            // slope for species
  //real b_Rhma;            // slope for species
  //real b_Osse;            // slope for species
  //real b_Raca;            // slope for species
  vector[K] r_micro;
  //real b_Anbo_inter;            // slope for interaction
  //real b_Rhma_inter;            // slope for interaction
  //real b_Osse_inter;            // slope for interaction
  //real b_Raca_inter;            // slope for interaction
  real b_micro;          // overall slope for microbiome predictor
  vector[N] x;          // unknown true value
  real<lower=0> sigma;  // sample-level variation
  real<lower=0> ksigma; // random effect variation
}
model {
  //b_Anbo ~ normal(0, 10);
  //b_Rhma ~ normal(0, 10);
  //b_Osse ~ normal(0, 10);
  //b_Raca ~ normal(0, 10);
  b_micro ~ normal(0,10);
  sigma ~ cauchy(0, 5);
  ksigma ~ cauchy(0,5);
  for ( k in 1:K ) 
      r_micro[k] ~ normal(b_micro, ksigma);
  
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ normal(x_species*b_species + x_species*r_micro + b_micro * x, sigma);
  
}
