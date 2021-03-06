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

//Test
vector[n] test_dose_mg_twice_daily;
vector[n] test_hrs_post_dose;


vector[n] test_age;
vector[n] test_sex;
vector[n] test_weight_kg;
vector[n] test_creatinine_micromol_l;

//Drug
vector[n] test_amiodarone_mg_day;
vector[n] test_diltiazem_mg_day;
}
transformed data{
  vector[n] yobs = yobs_ng_ml * 0.001;
  vector[n] scaled_creatinine = (creatinine_micromol_l - mean(creatinine_micromol_l))/sd(creatinine_micromol_l);
  vector[n] scaled_weight = (weight_kg - mean(weight_kg))/sd(weight_kg);
  vector[n] scaled_age = (age - mean(age))/sd(age);
  
  vector[n] scaled_diltiazem = diltiazem_mg_day/max(diltiazem_mg_day);
  vector[n] scaled_amiodarone = amiodarone_mg_day/max(amiodarone_mg_day);
  
  //Test
  vector[n] scaled_test_creatinine = (test_creatinine_micromol_l - mean(creatinine_micromol_l))/sd(creatinine_micromol_l);
  vector[n] scaled_test_weight = (test_weight_kg - mean(weight_kg))/sd(weight_kg);
  vector[n] scaled_test_age = (test_age - mean(age))/sd(age);
  
  vector[n] scaled_test_diltiazem = test_diltiazem_mg_day/max(diltiazem_mg_day);
  vector[n] scaled_test_amiodarone = test_amiodarone_mg_day/max(amiodarone_mg_day);
  
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
  
  real<lower=0> tmax;
  
  real mu_a;
  real b_sex_a;
  real<lower=0> sigma;
  
}
transformed parameters{
 vector<lower=0>[n] V = exp(mu_V + b_sex_V*sex);
 vector<lower=0>[n] ke = exp(mu_ke + b_creatinine_ke*scaled_creatinine + b_weight_ke*scaled_weight + b_age_ke*scaled_age + b_sex_ke*sex);
 vector<lower=0, upper=1>[n] F = inv_logit(mu_F + b_amiodarone_F*scaled_amiodarone + b_diltiazem_F*scaled_diltiazem);
 vector<lower=0, upper=1>[n] a = inv_logit(mu_a + b_sex_a*sex);
 vector<lower=0>[n] ka= -log(a)./(tmax*(1-a));
 vector<lower=0>[n] C0 = rep_vector(0.0, n);
 vector<lower=0>[n] C;
 
 for(i in 1:5){
   C0 = C0 + concentration(rep_vector(12*i, n), dose_mg_twice_daily, F, V, ka, ke);
 }
 
  C = C0 + concentration(hrs_post_dose, dose_mg_twice_daily, F, V, ka, ke);
}
model{
  mu_V ~ normal(log(21), .5);
  b_sex_V ~ normal(0, 0.075);
  
  mu_ke ~ normal(-1.25, .25);
  b_creatinine_ke ~ normal(0, 0.25);
  b_weight_ke ~ normal(0, 0.25);
  b_age_ke ~ normal(0, 0.25);
  b_sex_ke ~ normal(0, 0.25);
  
  tmax ~ normal(3.5, 0.5
  );
  
  mu_a ~ normal(-1, 0.25);
  b_sex_a ~ normal(0, 0.25);
  
  mu_F ~ normal(0, 0.125);
  b_amiodarone_F ~ double_exponential(0, tau_F);
  b_diltiazem_F ~ double_exponential(0, tau_F);
  tau_F ~ normal(0,1);
  
  sigma ~ lognormal(log(0.1), 0.2);
  
  yobs ~ lognormal(log(C), sigma);  
  
  
}
generated quantities{
 vector[n] yppc = 1000*to_vector(lognormal_rng(log(C), sigma));
  
  vector[n] test_V = exp(mu_V + b_sex_V*test_sex);
  vector<lower=0>[n] test_ke = exp(mu_ke + b_creatinine_ke*scaled_test_creatinine + b_weight_ke*scaled_test_weight + b_age_ke*scaled_test_age + b_sex_ke*test_sex);
  vector<lower=0, upper=1>[n] test_a = inv_logit(mu_a + b_sex_a*test_sex);
  vector<lower=0>[n] test_ka=  -log(test_a)./(tmax*(1-test_a));
  vector<lower=0, upper=1>[n] test_F = inv_logit(mu_F + b_amiodarone_F*scaled_test_amiodarone + b_diltiazem_F*scaled_test_diltiazem);
  vector[n] test_C0 = rep_vector(0.0, n);
  vector[n] test_C;
  vector[n] test_yppc;

  for(i in 1:5){
    test_C0 = test_C0 + concentration(rep_vector(12*i, n), test_dose_mg_twice_daily, test_F, test_V, test_ka, test_ke);
  }

  test_C = test_C0 + concentration(test_hrs_post_dose, test_dose_mg_twice_daily, test_F, test_V, test_ka, test_ke);

  test_yppc = 1000*to_vector(lognormal_rng(log(test_C), sigma));


}