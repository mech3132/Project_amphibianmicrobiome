//Stan normal model with uncertainty NO RAPI
data {
  int<lower=0> N;       // number of cases
  real<lower=0> y[N];          // outcome (variate)
  //vector[N] y;          // outcome (variate)
  vector[N] ysd;        //outcome uncertainty
  vector[N] x_meas;     // measurement of x
  vector[N] tau;       // sd of measurment of x
  vector[N] x_Anbo;       // species
  vector[N] x_Rhma;       // species
  vector[N] x_Osse;       // species
  vector[N] x_Raca;       // species
  //vector[N] x_Rapi;       // species
}
parameters {
  real b_Anbo;            // slope for species
  real b_Rhma;            // slope for species
  real b_Osse;            // slope for species
  real b_Raca;            // slope for species
  //real b_Rapi;            // slope for species
  real b_micro;          // slope for microbiome predictor
  vector[N] x;          // unknown true value
  //real<lower=0> sigma;  // sample-level variation
  //vector[N] y2;
}
model {
  b_Anbo ~ normal(0, 10);
  b_Rhma ~ normal(0, 10);
  b_Osse ~ normal(0, 10);
  b_Raca ~ normal(0, 10);
  //b_Rapi ~ normal(0, 10);
  b_micro ~ normal(0,10);
  //sigma ~ cauchy(0, 5);
  x ~ normal(0, 10);  // prior
  x_meas ~ normal(x, tau);    // measurement model
  y ~ normal(b_Anbo * x_Anbo + b_Rhma * x_Rhma + b_Osse * x_Osse + b_Raca * x_Raca + b_micro * x, ysd);
  //y ~ normal(y2, ysd);
  
}
