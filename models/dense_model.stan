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
}
transformed data{
  vector[r_n_subjectids] r_scaled_weight = (r_weight - mean(r_weight))/sd(r_weight);
  vector[r_n_subjectids] r_scaled_age = (r_age - mean(r_age))/sd(r_age);
  vector[r_n_subjectids] r_scaled_creatinine = (r_creatinine - mean(r_creatinine))/sd(r_creatinine);
  matrix[r_n_subjectids, 4] X = [r_is_male', r_scaled_weight', r_scaled_creatinine', r_scaled_age']';
  vector[r_n] r_yobs_scaled = r_yobs/1000;

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
  
  b_a_creatinine ~ normal(0, 0.25);
  
  r_sigma ~ lognormal(log(0.1), 0.2);
  r_yobs_scaled ~ lognormal(log(r_C), r_sigma);
}
generated quantities{
  real r_yppc[r_n] = lognormal_rng(log(r_C), r_sigma);
}

