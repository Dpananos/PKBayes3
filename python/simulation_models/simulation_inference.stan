functions{
  vector concentration(vector t, vector D, vector F, vector Cl, vector ka, vector ke){
    return D .* F .* ka .* ke ./ (Cl .* (ke - ka)) .* (exp(-ka .* t) - exp(-ke .* t));
  }
}
data{
  //Rommel's data
  int dense_n; //Total number of observations
  int dense_subjectids[dense_n]; //Subject idendification number as an integer.  Mine go from 1 - 36
  int dense_n_subjectids; //How many unique subjects do I have?
  vector[dense_n] dense_time; //time at which subjects were observed?  Length N
  vector[dense_n] dense_yobs; //Observed concentraitons
  
  vector[dense_n_subjectids] dense_is_male;
  vector[dense_n_subjectids] dense_weight;
  vector[dense_n_subjectids] dense_creatinine;
  vector[dense_n_subjectids] dense_age;
  vector[dense_n_subjectids] dense_dose;
  
  //Ute's data
  int sparse_n; //Total number of observations
  vector[sparse_n] sparse_time; //time at which subjects were observed?  Length N
  vector[sparse_n] sparse_yobs; //Observed concentraitons
  
  vector[sparse_n] sparse_is_male;
  vector[sparse_n] sparse_weight;
  vector[sparse_n] sparse_creatinine;
  vector[sparse_n] sparse_age;
  vector[sparse_n] sparse_dose;
  vector[sparse_n] sparse_amio;

}
parameters{
  
  real  mu_cl;
  real<lower=0> s_cl;                                                   
  vector[dense_n_subjectids] z_cl;
  
  real mu_tmax;
  real<lower=0> s_t;
  vector[dense_n_subjectids] z_t;
  
  real<lower=0> dense_sigma;
  real<lower=0> sparse_sigma;
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[dense_n_subjectids] z_alpha;
  
  real mu_F;
  
  real b_cl_age;
  real b_cl_weight;
  real b_cl_is_male;
  real b_cl_creatinine;
  
  real b_t_age;
  real b_t_weight;
  real b_t_is_male;
  real b_t_creatinine;

  real b_a_creatinine;
  
  real b_F_amio;

  real<lower=0> tau_F;
}
transformed parameters{

  vector<lower=0>[dense_n_subjectids] dense_Cl = exp(mu_cl + z_cl*s_cl + b_cl_age*dense_age + b_cl_weight*dense_weight + b_cl_is_male*dense_is_male + b_cl_creatinine*dense_creatinine );
  vector<lower=0>[dense_n_subjectids] dense_tmax = exp(mu_tmax + z_t*s_t + b_t_age*dense_age + b_t_weight*dense_weight + b_t_is_male*dense_is_male + b_t_creatinine*dense_creatinine);
  vector<lower=0, upper=1>[dense_n_subjectids] dense_alpha = inv_logit(mu_alpha + z_alpha*s_alpha + b_a_creatinine*dense_creatinine);
  vector<lower=0>[dense_n_subjectids] dense_ka = log(dense_alpha)./(dense_tmax .* (dense_alpha-1));
  vector<lower=0>[dense_n_subjectids] dense_ke = dense_alpha .* log(dense_alpha)./(dense_tmax .* (dense_alpha-1));
  real<lower=0, upper=1> dense_F = inv_logit(mu_F);  
  vector<lower=0>[dense_n] dense_C = concentration(dense_time, dense_dose[dense_subjectids], rep_vector(dense_F, dense_n), dense_Cl[dense_subjectids] ,  dense_ka[dense_subjectids], dense_ke[dense_subjectids]);

 vector<lower=0>[sparse_n] sparse_Cl = exp(mu_cl + b_cl_age*sparse_age + b_cl_weight*sparse_weight + b_cl_is_male*sparse_is_male + b_cl_creatinine*sparse_creatinine );
  vector<lower=0>[sparse_n] sparse_tmax = exp(mu_tmax + b_t_age*sparse_age + b_t_weight*sparse_weight + b_t_is_male*sparse_is_male + b_t_creatinine*sparse_creatinine );
  vector<lower=0, upper=1>[sparse_n] sparse_alpha = inv_logit(mu_alpha+ b_a_creatinine*sparse_creatinine);
  vector<lower=0>[sparse_n] sparse_ka = log(sparse_alpha)./(sparse_tmax .* (sparse_alpha-1));
  vector<lower=0>[sparse_n] sparse_ke = sparse_alpha .* sparse_ka;
  vector<lower=0, upper=1>[sparse_n] sparse_F = inv_logit(mu_F + b_F_amio*sparse_amio);

  vector<lower=0>[sparse_n] sparse_C = rep_vector(0.0, sparse_n);
  vector<lower=0>[sparse_n] sparse_C0 = rep_vector(0.0, sparse_n);

  for(i in 1:14){
    sparse_C0 = sparse_C0 + concentration(rep_vector(12*i, sparse_n), sparse_dose, sparse_F, sparse_Cl, sparse_ka, sparse_ke);
  }
  sparse_C = sparse_C0 + concentration(sparse_time, sparse_dose, sparse_F, sparse_Cl, sparse_ka, sparse_ke);

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
  
  b_cl_age ~ normal(0, 0.25);
  b_cl_weight ~ normal(0, 0.25);
  b_cl_is_male ~ normal(0, 0.25);
  b_cl_creatinine ~ normal(0, 0.25);
  
  b_t_age ~ normal(0, 0.25);
  b_t_weight ~ normal(0, 0.25);
  b_t_is_male ~ normal(0, 0.25);
  b_t_creatinine ~ normal(0, 0.25);
  
  b_a_creatinine ~ normal(0, 0.25);
  
  b_F_amio ~ double_exponential(0, tau_F);
  tau_F ~ normal(0, 0.25);
  
  dense_sigma ~ lognormal(log(0.1), 0.2);
  sparse_sigma ~ lognormal(log(0.1), 0.2);
  dense_yobs ~ lognormal(log(dense_C), dense_sigma);
  sparse_yobs ~ lognormal(log(sparse_C), sparse_sigma);
}
generated quantities{
  real dense_yppc[dense_n] = lognormal_rng(log(dense_C), dense_sigma);
  real sparse_yppc[sparse_n] = lognormal_rng(log(sparse_C), sparse_sigma);
}

