data {
  int<lower=0> N; // number of observations
  vector[N] y; // vector of observations
  vector<lower=0, upper=1>[N] x1; // (1 - g/2) * (1 - t)
  vector<lower=0, upper=1>[N] x2; // (g/2) * (1 - t)
  vector<lower=0, upper=1>[N] x3; // (1 - g/2) * t
  vector<lower=0, upper=1>[N] x4; // (g/2) * t
  int<lower=0, upper=1> m1; // {0,1} selecting coefficient of genotype
  int<lower=0, upper=1> m2; // {0,1} selecting coefficient of treatment
  int<lower=0, upper=1> m3; // {0,1} selecting coefficient of interaction
  real<lower=0> phi1; // scaled prior standard deviation of effects
  real<lower=0> phi2; // scaled prior standard deviation of effects
  real<lower=0> phi3; // scaled prior standard deviation of effects
}
parameters {
  real b0; // intercept
  real b1; // coefficient of genotype
  real b2; // coefficient of treatment
  real b3; // coefficient of interaction
  real<lower=0> sigma_sq_inv; // precision of error
}
transformed parameters {
  real<lower=0> sigma;
  sigma = sigma_sq_inv^(-0.5); // standard deviation of error
}
model {
  vector[N] r_mu; // note that constraint cannot be put in the model block
  r_mu = x1 + exp(2*m1*b1) * x2 + exp(m2*b2) * x3 + exp(2*m1*b1 + m2*b2 + 2*m3*b3) * x4;
  target += gamma_lpdf(sigma_sq_inv | 0.001, 0.001); // prior for precision of error
  target += normal_lpdf(b0 | 0, sigma*sqrt(1000)); // prior for intercept
  target += normal_lpdf(b1 | 0, sigma*phi1); // prior for coefficient of genotype; untruncated
  target += normal_lpdf(b2 | 0, sigma*phi2); // prior for coefficient of treatment; untruncated
  target += normal_lpdf(b3 | 0, sigma*phi3); // prior for coefficient of interaction; untruncated
  target += normal_lpdf(y | b0 + log(r_mu), sigma); // likelihood model
}
