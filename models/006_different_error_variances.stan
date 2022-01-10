functions{
  vector concentration(vector t, vector D, vector F, vector V, vector ka, vector ke){
    return D .* F .* ka ./ (V .* (ke - ka)) .* (exp(-ka .* t) - exp(-ke .* t));
  }
}
data{
  //Basic
  int n;
  int p_v;
  int p_ke;
  int p_F;
  int p_a;
  
  matrix[n, p_v] X_v;
  matrix[n, p_ke] X_ke;
  matrix[n, p_F] X_F;
  matrix[n, p_a] X_a;
  
  vector[n] dose_mg_twice_daily;
  vector[n] hrs_post_dose;
  vector[n] yobs_ng_ml;
  
  
}
transformed data{
  vector[n] yobs = yobs_ng_ml * 0.001;
  vector[n] dose_size = (dose_mg_twice_daily-2.5)/2.5;
}
parameters{
  vector[p_v] beta_v;
  vector[p_ke] beta_ke;
  vector[p_F] beta_F;
  vector[p_a] beta_a;
  real<lower=0> tmax;
  real sigma_0;
  real sigma_1;
  real<lower=0>tau_v;
  real<lower=0>tau_ke;
  real<lower=0>tau_F;
  real<lower=0>tau_a;
}
transformed parameters{
  vector<lower=0>[n] V = exp(X_v * beta_v);
  vector<lower=0>[n] ke = exp(X_ke * beta_ke);
  vector<lower=0, upper=1>[n] f = inv_logit(X_F * beta_F);
  vector<lower=0, upper=1>[n] a = inv_logit(X_a * beta_a);
  vector<lower=0>[n] ka = -log(a)./(tmax*(1-a));
  vector<lower=0>[n] C0 = rep_vector(0.0, n);
  vector<lower=0>[n] C;
  vector<lower=0>[n] sigma = exp(sigma_0 + dose_size*sigma_1);
  
  for(i in 1:10){
    C0 = C0 + concentration(rep_vector(12*i, n), dose_mg_twice_daily, f, V, ka, ke);
  }
  
  C = C0 + concentration(hrs_post_dose, dose_mg_twice_daily, f, V, ka, ke);
}
model{
  beta_v[1] ~ normal(log(21.0), 1.0);
  beta_v[2:p_v] ~ double_exponential(0, tau_v);
  
  beta_ke[1] ~ normal(-1.25, 0.25);
  beta_ke[2:p_ke] ~ double_exponential(0, tau_ke);
  
  beta_a[1] ~ normal(-1, 0.25);
  beta_a[2:p_a] ~ double_exponential(0, tau_a);
  
  beta_F[1] ~ lognormal(log(0.1), 0.2);
  beta_F[2:p_F] ~ double_exponential(0, tau_F);
  
  tmax ~ normal(3.5, 0.5);
  
  tau_v ~ normal(0, 1);
  tau_ke ~ normal(0, 1);
  tau_F ~ normal(0, 1);
  tau_a ~ normal(0, 1);
  sigma_0 ~  normal(log(0.1), 0.05);
  sigma_1 ~ normal(0, 1);
  yobs ~ lognormal(log(C), sigma);
}
generated quantities{
  vector[n] yppc = 1000*to_vector(lognormal_rng(log(C), sigma));
}