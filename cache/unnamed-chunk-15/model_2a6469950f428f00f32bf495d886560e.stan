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
}
parameters{
  vector[p_v] beta_v;
  vector[p_ke] beta_ke;
  vector[p_F] beta_F;
  vector[p_a] beta_a;
  real<lower=0> sigma;
}
transformed parameters{
 vector[n] V = exp(X_v * beta_v);
 vector[n] ke = exp(X_ke * beta_ke);
 vector[n] f = inv_logit(X_F * beta_F);
 vector[n] a = inv_logit(X_a * beta_a);
 vector[n] ka = ke ./ a;
 vector[n] C = concentration(hrs_post_dose, dose_mg_twice_daily, f, V, ka, ke);
}
model{
  beta_v ~ normal(log(21.0), 1.0);
  beta_ke ~ normal(0, 1);
  beta_a ~ normal(0, 0.125);
  beta_F ~ normal(0, 0.125);
  sigma ~ cauchy(0, 1);
  yobs ~ lognormal(log(C), sigma);
}
generated quantities{
  vector[n] yppc = 1000*to_vector(lognormal_rng(log(C), sigma));
}
