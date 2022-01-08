functions{
  vector concentration(vector t, vector D, vector F, vector Cl, vector ka, vector ke){
    return D .* F .* ka .* ke ./ (Cl .* (ke - ka)) .* (exp(-ka .* t) - exp(-ke .* t));
  }
}
data{
  //Rommel's data
  int r_n; //Total number of observations
  int r_subjectids[r_n]; //Subject idendification number as an integer.  Mine go from 1 - 36
  int r_n_subjectids; //How many unique subjects do I have?
  vector[r_n] r_time; //time at which subjects were observed?  Length N
  vector[r_n] r_yobs; //Observed concentraitons
  
  vector[r_n_subjectids] r_is_male;
  vector[r_n_subjectids] r_weight;
  vector[r_n_subjectids] r_creatinine;
  vector[r_n_subjectids] r_age;
  vector[r_n_subjectids] r_D;
  
  //Ute's data
  int u_n; //Total number of observations
  vector[u_n] u_time; //time at which subjects were observed?  Length N
  vector[u_n] u_yobs; //Observed concentraitons
  
  vector[u_n] u_is_male;
  vector[u_n] u_weight;
  vector[u_n] u_creatinine;
  vector[u_n] u_age;
  vector[u_n] u_D;
  vector[u_n] u_amiodarone;
  vector[u_n] u_diltiazem;

}
transformed data{
  vector[r_n_subjectids] r_scaled_weight = (r_weight - mean(r_weight))/sd(r_weight);
  vector[r_n_subjectids] r_scaled_age = (r_age - mean(r_age))/sd(r_age);
  vector[r_n_subjectids] r_scaled_creatinine = (r_creatinine - mean(r_creatinine))/sd(r_creatinine);
  matrix[r_n_subjectids, 4] X = [r_is_male', r_scaled_weight', r_scaled_creatinine', r_scaled_age']';
  vector[r_n] r_yobs_scaled = r_yobs/1000;
  
  
  vector[u_n] u_scaled_weight = (u_weight - mean(r_weight))/sd(r_weight);
  vector[u_n] u_scaled_age = (u_age - mean(r_age))/sd(r_age);
  vector[u_n] u_scaled_creatinine = (u_creatinine - mean(r_creatinine))/sd(r_creatinine);
  vector[u_n] u_yobs_scaled = u_yobs/1000;
  matrix[u_n, 4] u_X = [u_is_male', u_scaled_weight', u_scaled_creatinine', u_scaled_age']';
  
  vector[u_n] u_amiodarone_scaled = u_amiodarone/max(u_amiodarone);
  vector[u_n] u_diltiazem_scaled = u_diltiazem/max(u_diltiazem);

}
parameters{
  
  real  mu_cl;
  real<lower=0> s_cl;                                                   
  vector[r_n_subjectids] z_cl;
  
  real mu_tmax;
  real<lower=0> s_t;
  vector[r_n_subjectids] z_t;
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[r_n_subjectids] delays;
  
  real<lower=0> sigma;
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[r_n_subjectids] z_alpha;
  
  real mu_F;
  
  vector[4] beta_cl;
  vector[4] beta_t;
  vector[4] beta_a;
  real beta_amio;
  real beta_dil;
  real<lower=0> tau_F;
}
transformed parameters{
  //Rommel's parameters
  vector<lower=0>[r_n_subjectids] r_Cl = exp(mu_cl + z_cl*s_cl + X*beta_cl);
  vector<lower=0>[r_n_subjectids] r_tmax = exp(mu_tmax + z_t*s_t + X*beta_t);
  vector<lower=0, upper=1>[r_n_subjectids] r_alpha = inv_logit(mu_alpha + z_alpha*s_alpha + X*beta_a);
  vector<lower=0>[r_n_subjectids] r_ka = log(r_alpha)./(r_tmax .* (r_alpha-1));
  vector<lower=0>[r_n_subjectids] r_ke = r_alpha .* log(r_alpha)./(r_tmax .* (r_alpha-1));
  vector<lower=0>[r_n] delayed_time = r_time - 0.5*delays[r_subjectids];
  real<lower=0, upper=1> r_F = inv_logit(mu_F);
  
  vector<lower=0>[r_n] r_C = concentration(delayed_time, r_D[r_subjectids], rep_vector(r_F, r_n), r_Cl[r_subjectids] ,  r_ka[r_subjectids], r_ke[r_subjectids]);

  //Ute's Parameters
  vector<lower=0>[u_n] u_Cl = exp(mu_cl + u_X*beta_cl);
  vector<lower=0>[u_n] u_tmax = exp(mu_tmax + u_X*beta_t);
  vector<lower=0, upper=1>[u_n] u_alpha = inv_logit(mu_alpha + u_X*beta_a);
  vector<lower=0>[u_n] u_ka = log(u_alpha)./(u_tmax .* (u_alpha-1));
  vector<lower=0>[u_n] u_ke = u_alpha .* log(u_alpha)./(u_tmax .* (u_alpha-1));
  vector<lower=0, upper=1>[u_n] u_F = inv_logit(mu_F + beta_dil*u_diltiazem_scaled + beta_amio*u_amiodarone_scaled);

  vector<lower=0>[u_n] u_C = rep_vector(0.0, u_n);

  for(i in 1:10){
    u_C = u_C + concentration(rep_vector(12*i, u_n), u_D, u_F, u_Cl, u_ka, u_ke);
  }
  u_C = u_C + concentration(u_time, u_D, u_F, u_Cl, u_ka, u_ke);

}
model{
  //See Byon et. al 2019
  mu_tmax ~ normal(log(3.3), 0.1);
  s_t ~ gamma(5, 100);
  z_t ~ normal(0,1);
  
  mu_cl ~ normal(log(3.3),0.15);
  s_cl ~ gamma(15,100);
  z_cl ~ normal(0,1);
  
  
  mu_alpha ~ normal(-0.25,0.5);
  s_alpha ~ gamma(10, 100);
  z_alpha ~ normal(0,1);
  
  mu_F ~ normal(0, 0.125);

  phi ~ beta(20,20);
  kappa ~ beta(20,20);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
  beta_cl ~ normal(0,0.25);
  beta_t ~ normal(0, 0.25);
  beta_a ~ normal(0, 0.25);
  beta_amio ~ double_exponential(0, tau_F);
  beta_dil ~ double_exponential(0, tau_F);
  tau_F ~ normal(0, 0.25);
  
  sigma ~ lognormal(log(0.1), 0.2);
  r_yobs_scaled ~ lognormal(log(r_C), sigma);
  u_yobs_scaled ~ lognormal(log(u_C), sigma);
}
generated quantities{
  real r_yppc[r_n] = lognormal_rng(log(r_C), sigma);
  real u_yppc[u_n] = lognormal_rng(log(u_C), sigma);
}
