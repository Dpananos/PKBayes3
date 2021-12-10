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
}
parameters{
  real<lower=0> V;
  real<lower=0> ke;
  real<lower=0, upper=1> F;
  real<lower=0, upper=1> a;
  real<lower=0> sigma;
}
transformed parameters{
 real<lower=0> ka= ke/a;
 vector[n] C = concentration(hrs_post_dose, dose_mg_twice_daily, rep_vector(F, n), rep_vector(V, n), rep_vector(ka, n), rep_vector(ke, n));
}
model{
  V ~ normal(21.0, 1.0);
  ke ~ normal(0, 1);
  a ~ beta(2, 2);
  F ~ normal(0.5, 0.05);
  sigma ~ cauchy(0, 1);
  yobs ~ lognormal(log(C), sigma);
}
generated quantities{
  vector[n] yppc = 1000*to_vector(lognormal_rng(log(C), sigma));
  vector[n] Cpred = 1000*C;
}