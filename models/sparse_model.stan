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
  
  real mu_tmax;
  
  
  real<lower=0> u_sigma;
  
  real mu_alpha;
  
  real mu_F;
  
  real b_cl_age;
  real b_cl_weight;
  real b_cl_is_male;
  real b_cl_creatinine;
  
  real b_t_age;
  real b_t_weight;
  real b_t_is_male;
  real b_t_creatinine;
  
  // real b_a_age;
  // real b_a_weight;
  // real b_a_is_male;
  real b_a_creatinine;
  
  real b_F_amio;
  // real beta_dil;
  real<lower=0> tau_F;
}
transformed parameters{
  //Rommel's parameters

  //Ute's Parameters
  vector<lower=0>[u_n] u_Cl = exp(mu_cl + b_cl_age*u_scaled_age + b_cl_weight*u_scaled_weight + b_cl_is_male*u_is_male + b_cl_creatinine*u_scaled_creatinine );
  vector<lower=0>[u_n] u_tmax = exp(mu_tmax + b_t_age*u_scaled_age + b_t_weight*u_scaled_weight + b_t_is_male*u_is_male + b_t_creatinine*u_scaled_creatinine );
  vector<lower=0, upper=1>[u_n] u_alpha = inv_logit(mu_alpha+ b_a_creatinine*u_scaled_creatinine);
  vector<lower=0>[u_n] u_ka = log(u_alpha)./(u_tmax .* (u_alpha-1));
  vector<lower=0>[u_n] u_ke = u_alpha .* u_ka;
  vector<lower=0, upper=1>[u_n] u_F = inv_logit(mu_F + b_F_amio*u_amiodarone_scaled);
  
  vector<lower=0>[u_n] u_C = rep_vector(0.0, u_n);
  vector<lower=0>[u_n] u_C0 = rep_vector(0.0, u_n);
  
  for(i in 1:14){
    u_C0 = u_C0 + concentration(rep_vector(12*i, u_n), u_D, u_F, u_Cl, u_ka, u_ke);
  }
  u_C = u_C0 + concentration(u_time, u_D, u_F, u_Cl, u_ka, u_ke);
  
  
  
}
model{
  //See Byon et. al 2019
  mu_tmax ~ normal(log(3.3), 0.1);
  
  mu_cl ~ normal(log(3.3),0.15);
  
  
  mu_alpha ~ normal(-0.25,0.5);
  
  mu_F ~ normal(0, 0.025);

  
  b_cl_age ~ normal(0, 0.25);
  b_cl_weight ~ normal(0, 0.25);
  b_cl_is_male ~ normal(0, 0.25);
  b_cl_creatinine ~ normal(0, 0.25);
  
  b_t_age ~ normal(0, 0.25);
  b_t_weight ~ normal(0, 0.25);
  b_t_is_male ~ normal(0, 0.25);
  b_t_creatinine ~ normal(0, 0.25);
  
  // b_a_age ~ normal(0, 0.25);
  // b_a_weight ~ normal(0, 0.25);
  // b_a_is_male ~ normal(0, 0.25);
  b_a_creatinine ~ normal(0, 0.25);
  
  b_F_amio ~ double_exponential(0, tau_F);
  // beta_dil ~ double_exponential(0, tau_F);
  tau_F ~ normal(0, 0.25);
  
  u_sigma ~ lognormal(log(0.1), 0.2);
  u_yobs_scaled ~ lognormal(log(u_C), u_sigma);
}
generated quantities{
  real u_yppc[u_n] = lognormal_rng(log(u_C), u_sigma);
}

