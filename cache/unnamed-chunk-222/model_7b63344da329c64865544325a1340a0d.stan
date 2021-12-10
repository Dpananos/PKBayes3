functions{
  vector concentration(vector t, vector D, vector F, vector V, vector ka, vector ke){
   return D .* F .* ka ./ (V .* (ke - ka)) .* (exp(-ka .* t) - exp(-ke .* t));
  }
}
data{
//Basic
int n;
vector[n] dose_mg_twice_daily;
vector[n] hrs_post_dose;
vector[n] yobs_ng_ml;

//Clinical
vector[n] age;
vector[n] sex;
vector[n] weight_kg;
vector[n] creatinine_micromol_l;

//History
vector[n] history_arrhythmia_treatments_e_g_ablation_cardioversion;
vector[n] history_diabetes;
vector[n] history_hypertension;
vector[n] history_ischemic_stroke;
vector[n] history_tia;
vector[n] history_hemorrhage;
vector[n] history_myocardial_infarction;
vector[n] history_pci;
vector[n] history_coronary_artery_bypass_graft;
vector[n] history_peripheral_vascular_disease;
vector[n] history_cancer;
vector[n] history_ckd;

//Drug
vector[n] amiodarone_mg_day;
vector[n] carbamazepine_mg_day;
vector[n] diltiazem_mg_day;
vector[n] phenobarbital_mg_day;
vector[n] phenytoin_mg_day;
vector[n] primidone_mg_day;
vector[n] rifampin_mg_day;
vector[n] verapamil_mg_day;
vector[n] ketoconazole_mg_day;
vector[n] fluconazole_mg_day;
}
transformed data{
  vector[n] yobs = yobs_ng_ml * 0.001;
  vector[n] scaled_creatinine = (creatinine_micromol_l - mean(creatinine_micromol_l))/sd(creatinine_micromol_l);
  vector[n] scaled_weight = (weight_kg - mean(weight_kg))/sd(weight_kg);
  vector[n] scaled_age = (age - mean(age))/sd(age);
  
  vector[n] scaled_diltiazem = diltiazem_mg_day/max(diltiazem_mg_day);
  vector[n] scaled_amiodarone = amiodarone_mg_day/max(amiodarone_mg_day);
}
parameters{
  
  real<lower=0> mu_V;
  real b_sex_V;
  
  real<lower=0> mu_ke;
  real b_creatinine_ke;
  real b_weight_ke;
  real b_age_ke;
  real b_sex_ke;
  
  
  real mu_F;
  real<lower=0> tau_F;
  real b_amiodarone_F;
  real b_diltiazem_F;
  
  real<lower=0, upper=1> a;
  real<lower=0> sigma;
}
transformed parameters{
 vector<lower=0>[n] V = mu_V*exp(b_sex_V*sex);
 vector<lower=0>[n] ke = mu_ke*exp(b_creatinine_ke*scaled_creatinine + b_weight_ke*scaled_weight + b_age_ke*scaled_age + b_sex_ke*sex);
 vector<lower=0, upper=1>[n] F = inv_logit(mu_F + b_amiodarone_F*scaled_amiodarone + b_diltiazem_F*scaled_diltiazem);
 vector<lower=0>[n] ka= ke/a;
 vector<lower=0>[n] C0 = rep_vector(0.0, n);
 vector[n] C;
 
 for(i in 1:5){
   C0 = C0 + concentration(rep_vector(12*i, n), dose_mg_twice_daily, F, V, ka, ke);
 }
 
  C = C0 + concentration(hrs_post_dose, dose_mg_twice_daily, F, V, ka, ke);
}
model{
  mu_V ~ normal(21.0, 1.0);
  b_sex_V ~ normal(0, 0.25);
  
  mu_ke ~ normal(0, 1);
  b_creatinine_ke ~ normal(0, 0.25);
  b_weight_ke ~ normal(0, 0.25);
  b_age_ke ~ normal(0, 0.25);
  b_sex_ke ~ normal(0, 0.25);
  
  a ~ beta(2, 2);
  
  mu_F ~ normal(0, 0.125);
  b_amiodarone_F ~ double_exponential(0, tau_F);
  b_diltiazem_F ~ double_exponential(0, tau_F);
  tau_F ~ normal(0,1);
  
  sigma ~ cauchy(0, 1);
  yobs ~ lognormal(log(C), sigma);
}
generated quantities{
  vector[n] yppc = 1000*to_vector(lognormal_rng(log(C), sigma));
}
