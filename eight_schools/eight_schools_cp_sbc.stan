functions {
  real positive_cauchy_rng(real mu, real sigma) {
    real p = cauchy_cdf(0, mu, sigma);
    real u = uniform_rng(p, 1);
    real y = mu + sigma*tan(pi()*(u - 0.5));
    return y;
  }
}
  
data {
  int<lower=0> J;
  real<lower=0> sigma[J];
}

transformed data {
  real mu_ = normal_rng(0,5);
  real tau_ = positive_cauchy_rng(0, 5);
  real theta_[J];
  real y[J];
  for (idx in 1:J) theta_[idx] = normal_rng(mu_, tau_);
  y = normal_rng(theta_, sigma);
}

parameters {
  real mu;
  real<lower=0> tau;
  real theta[J];
}

model {
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 5);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}

generated quantities {
  real y_[J] = y;
  vector[2+J] pars_;
  int ranks_[2+J];
  pars_[1] = mu_;
  pars_[2] = tau_;
  ranks_[1] = mu > mu_;
  ranks_[2] = tau > tau_;
  for (idx in 1:J) pars_[2+idx] = theta_[idx];
  for (idx in 1:J) ranks_[2+idx] = theta[idx] > theta_[idx];
}
