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
  
  int<lower=0> n_basis;
  matrix[r_n_subjectids, n_basis] r_age_spline;
  matrix[r_n_subjectids, n_basis] r_weight_spline;
  matrix[r_n_subjectids, n_basis] r_creat_spline;
  
  matrix[u_n, n_basis] u_age_spline;
  matrix[u_n, n_basis] u_weight_spline;
  matrix[u_n, n_basis] u_creat_spline;

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
  
  real<lower=0> r_sigma;
  real<lower=0> u_sigma;
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[r_n_subjectids] z_alpha;
  
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
  
  real beta_amio;
  // real beta_dil;
  real<lower=0> tau_F;
}
transformed parameters{
  //Rommel's parameters
  vector<lower=0>[r_n_subjectids] r_Cl = exp(mu_cl + z_cl*s_cl + b_cl_age*r_scaled_age + b_cl_weight*r_scaled_weight + b_cl_is_male*r_is_male + b_cl_creatinine*r_scaled_creatinine );
  vector<lower=0>[r_n_subjectids] r_tmax = exp(mu_tmax + z_t*s_t + b_t_age*r_scaled_age + b_t_weight*r_scaled_weight + b_t_is_male*r_is_male + b_t_creatinine*r_scaled_creatinine );
  vector<lower=0, upper=1>[r_n_subjectids] r_alpha = inv_logit(mu_alpha + z_alpha*s_alpha + b_a_creatinine*r_scaled_creatinine);
  vector<lower=0>[r_n_subjectids] r_ka = log(r_alpha)./(r_tmax .* (r_alpha-1));
  vector<lower=0>[r_n_subjectids] r_ke = r_alpha .* log(r_alpha)./(r_tmax .* (r_alpha-1));
  vector<lower=0>[r_n] delayed_time = r_time - 0.5*delays[r_subjectids];
  real<lower=0, upper=1> r_F = inv_logit(mu_F);
  
  vector<lower=0>[r_n] r_C = concentration(delayed_time, r_D[r_subjectids], rep_vector(r_F, r_n), r_Cl[r_subjectids] ,  r_ka[r_subjectids], r_ke[r_subjectids]);

  //Ute's Parameters
 vector<lower=0>[u_n] u_Cl = exp(mu_cl + b_cl_age*u_scaled_age + b_cl_weight*u_scaled_weight + b_cl_is_male*u_is_male + b_cl_creatinine*u_scaled_creatinine );
  vector<lower=0>[u_n] u_tmax = exp(mu_tmax + b_t_age*u_scaled_age + b_t_weight*u_scaled_weight + b_t_is_male*u_is_male + b_t_creatinine*u_scaled_creatinine );
  vector<lower=0, upper=1>[u_n] u_alpha = inv_logit(mu_alpha+ b_a_creatinine*u_scaled_creatinine);
  vector<lower=0>[u_n] u_ka = log(u_alpha)./(u_tmax .* (u_alpha-1));
  vector<lower=0>[u_n] u_ke = u_alpha .* u_ka;
  vector<lower=0, upper=1>[u_n] u_F = inv_logit(mu_F + beta_amio*u_amiodarone_scaled);

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
  s_t ~ gamma(5, 100);
  z_t ~ std_normal();
  
  mu_cl ~ normal(log(3.3),0.15);
  s_cl ~ gamma(15,100);
  z_cl ~ std_normal();
  
  
  mu_alpha ~ normal(-0.25,0.5);
  s_alpha ~ gamma(10, 100);
  z_alpha ~ normal(0,1);
  
  mu_F ~ normal(0, 0.025);

  phi ~ beta(20,20);
  kappa ~ beta(20,20);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
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
  
  beta_amio ~ double_exponential(0, tau_F);
  // beta_dil ~ double_exponential(0, tau_F);
  tau_F ~ normal(0, 0.25);
  
  r_sigma ~ lognormal(log(0.1), 0.2);
  u_sigma ~ lognormal(log(0.1), 0.2);
  r_yobs_scaled ~ lognormal(log(r_C), r_sigma);
  u_yobs_scaled ~ lognormal(log(u_C), u_sigma);
}
generated quantities{
  real r_yppc[r_n] = lognormal_rng(log(r_C), r_sigma);
  real u_yppc[u_n] = lognormal_rng(log(u_C), u_sigma);
}

