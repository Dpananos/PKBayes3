

data{
  int n;
  int n_subjectids;
  int subjectids[n];
  vector[n] time;
  vector[n] age;
  vector[n] weight;
  vector[n] creatinine;
  int is_male[n];
  vector[n] dose;

  real amio_effect;
}
generated quantities {
  // From previous paper
   real mu_cl = normal_rng(0.5, 0.04);
   real mu_t = normal_rng(0.93, 0.05);
   real mu_a = normal_rng(-1.35, 0.14);

   real s_cl = gamma_rng(69, 338);
   real s_t = gamma_rng(75, 350);
   real s_a = gamma_rng(10, 102);

   real z_cl = normal_rng(rep_vector(0, n_subjectids), 1);
   real z_t = normal_rng(rep_vector(0, n_subjectids), 1);
   real z_a = normal_rng(rep_vector(0, n_subjectids), 1);

   real b_cl_age = normal_rng(-0.15, 0.01);
   real b_cl_weight = normal_rng(0.09, 0.02);
   real b_cl_is_male = normal_rng(0.23, 0.04);
   real b_cl_creatinine = normal_rng(-0.07, 0.01);

   
}