functions{
vector system(vector y, vector theta, real[] x_r, int[] x_i ){
  real D = x_r[1];
  vector[3] z;
  
  
  
  
  z[1] = theta[1] - sqrt(exp(y[3])/exp(y[2]));
  z[2] = theta[2] - D*exp(y[1] - 2*sqrt(exp(y[2])*exp(y[3])));
  z[3] = theta[3] - sqrt(exp(y[2])/exp(y[3])) / (2*exp(y[1]) * bessel_second_kind(1, 2*sqrt(exp(y[2])*exp(y[3]))));
  
  return z;
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
  vector[3] b_guess = to_vector({1.5, 1.1, 1.9});
  vector[n] yobs = yobs_ng_ml * 0.001;
  int x_i[1] = {0};
  real x_r[1];
}
parameters{
 real<lower=0> tmax;
 real<lower=0> cmax;
 real<lower=0> cl;
 real<lower=0> phi;
}
transformed parameters{
  vector[3] theta = to_vector({tmax, cmax, cl});
  vector[3] b;
  vector[n] E_y;
  
  for(i in 1:n){
    
    b = algebra_solver(system, b_guess, theta, {dose_mg_twice_daily[i]}, x_i, 1e-10, 1e-6, 1e6); 
    E_y[i] = dose_mg_twice_daily[i] * exp(b[1] -exp(b[2])*hrs_post_dose[i] -exp(b[3])/hrs_post_dose[i]);
  }
}
model{
  tmax ~ normal(3.5, 0.25);
  cmax ~ normal(0.2, 0.05);
  cl ~ normal(3.3, 1.0);
  phi ~ cauchy(0, 1);
  for(i in 1:n){
    target += gamma_lpdf(yobs[i] | 1/phi, 1/(E_y[i]*phi));
  }
}
