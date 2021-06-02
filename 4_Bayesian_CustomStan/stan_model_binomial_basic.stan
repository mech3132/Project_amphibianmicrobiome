//Simple Stan bernoulli model
data {
  int<lower=0> N;       // number of cases
  int<lower=0, upper=1> y[N];          // outcome (variate)
  vector[N] x_meas;     // measurement of x
  vector[N] x_Anbo;       // species
  vector[N] x_Rhma;       // species
  vector[N] x_Osse;       // species
  vector[N] x_Raca;       // species
  vector[N] x_Rapi;       // species
  //real<lower=0> sdMean[N];  // measurement noise
  //real<lower=0> seMean[N]; // se of true mean
}
parameters {
  real alpha;           // intercept
  real b_Anbo;            // slope for species
  real b_Rhma;            // slope for species
  real b_Osse;            // slope for species
  real b_Raca;            // slope for species
  real b_Rapi;            // slope for species
  real b_micro;          // slope for microbiome predictor
  //real<lower=0> sigma;  // outcome noise
  //real x[N];          // unknown true value
  //real mu_x;          // prior location
  //real sigma_x;       // prior scale
}
model {
  y ~ bernoulli_logit(alpha + b_Anbo * x_Anbo + b_Rhma * x_Rhma + b_Osse * x_Osse + b_Raca * x_Raca+ b_Rapi * x_Rapi + b_micro * x_meas);
  alpha ~ normal(0, 10);
  b_Anbo ~ normal(0, 10);
  b_Rhma ~ normal(0, 10);
  b_Osse ~ normal(0, 10);
  b_Raca ~ normal(0, 10);
  b_Rapi ~ normal(0, 10);
  b_micro ~ normal(0,10);
  //sigma ~ cauchy(0, 5);
  //x ~ normal(mu_x, seMean);  // prior
  //x_meas ~ normal(x, sdMean);    // measurement model
  //y ~ normal(alpha + beta * x, sigma);
}
