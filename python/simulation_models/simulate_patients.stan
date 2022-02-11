functions{
  real pkprofile(real t, real D, real f, real cl, real ka, real ke){
    return D * f * ka * ke / (cl * (ke - ka)) * (exp(-ka * t) - exp(-ke * t));
  }
}

data{
  int n;
  int n_subjectids;
  int subjectids[n];
  vector[n] time;
  vector[n] age;
  vector[n] weight;
  vector[n] creatinine;
  vector[n] amio;
  vector[n] is_male;
  vector[n] dose;

  real amio_effect;
}
generated quantities {
  // From previous paper
   real mu_cl = normal_rng(0.5, 0.04);
   real mu_t = normal_rng(0.93, 0.05);
   real mu_a = normal_rng(-1.35, 0.14);
   real mu_f = normal_rng(0, 0.025);

   real s_cl = gamma_rng(69, 338);
   real s_t = gamma_rng(75, 350);
   real s_a = gamma_rng(10, 102);

   vector[n_subjectids] z_cl = to_vector(normal_rng(rep_vector(0, n_subjectids), 1));
   vector[n_subjectids] z_t = to_vector(normal_rng(rep_vector(0, n_subjectids), 1));
   vector[n_subjectids] z_a = to_vector(normal_rng(rep_vector(0, n_subjectids), 1));

   real b_cl_age = normal_rng(-0.15, 0.01);
   real b_cl_weight = normal_rng(0.09, 0.02);
   real b_cl_is_male = normal_rng(0.23, 0.04);
   real b_cl_creatinine = normal_rng(-0.07, 0.01);

   real b_t_age = normal_rng(0.03, 0.03);
   real b_t_weight = normal_rng(0.04, 0.04);
   real b_t_is_male = normal_rng(0.02, 0.06);
   real b_t_creatinine = normal_rng(-0.05, 0.05);

   real b_a_creatinine = normal_rng(-0.3, 0.1);

   vector[n] cl = exp(mu_cl + s_cl*z_cl[subjectids] + b_cl_age*age + b_cl_weight*weight + b_cl_is_male*is_male + b_cl_creatinine*creatinine);
   vector[n] tmax = exp(mu_t + s_t*z_t[subjectids] + b_t_age*age + b_t_weight*weight + b_t_is_male*is_male + b_t_creatinine*creatinine);
   vector[n] a = inv_logit(mu_a + s_a*z_a[subjectids] + b_a_creatinine*creatinine);
   vector[n] ka = log(a) ./ (tmax .* (a-1));
   vector[n] ke = a .* ka;
   vector[n] f = inv_logit(mu_f + amio_effect * amio);
   real c0 = 0.0;
   vector[n] concentration = rep_vector(0.0, n);
   real sigma = lognormal_rng(-1.75, 0.05);
   real observed_concentration[n];

   for(ii in 1:n){

     if(time[ii]<13){
       concentration[ii] = pkprofile(time[ii], dose[ii], f[ii], cl[ii], ka[ii], ke[ii]);
     }
     else{
      // Determine the initial condition
      // Simulate the concentration after taking 14 doses (~7 days into the dose regiment)
       c0 = 0.0;
       for(jj in 1:14){
        c0 = c0 + pkprofile(12*jj, dose[ii], f[ii], cl[ii], ka[ii], ke[ii]);
       }
       
       concentration[ii] = c0 + pkprofile(time[ii], dose[ii], f[ii], cl[ii], ka[ii], ke[ii]);
     }
   }

   observed_concentration = lognormal_rng(log(concentration), sigma);
}