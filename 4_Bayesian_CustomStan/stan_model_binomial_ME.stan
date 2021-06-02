//Stan bernoulli model with uncertainty
data {
  int<lower=0> N;       // number of cases
  int<lower=0, upper=1> y[N];          // outcome (variate)
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  vector[N] x_Anbo;       // species
  vector[N] x_Rhma;       // species
  vector[N] x_Osse;       // species
  vector[N] x_Raca;       // species
  vector[N] x_Rapi;       // species
}
parameters {
  real b_Anbo;            // slope for species
  real b_Rhma;            // slope for species
  real b_Osse;            // slope for species
  real b_Raca;            // slope for species
  real b_Rapi;            // slope for species
  real b_micro;          // slope for microbiome predictor
  vector[N] x;          // unknown true value
}
model {
  b_Anbo ~ normal(0, 10);
  b_Rhma ~ normal(0, 10);
  b_Osse ~ normal(0, 10);
  b_Raca ~ normal(0, 10);
  b_Rapi ~ normal(0, 10);
  b_micro ~ normal(0,10);
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ bernoulli_logit(b_Anbo * x_Anbo + b_Rhma * x_Rhma + b_Osse * x_Osse + b_Raca * x_Raca+ b_Rapi * x_Rapi + b_micro * x);
  
}
