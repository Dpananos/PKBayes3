library(tidyverse)
library(posterior)
library(tidybayes)

prepare_bayesian_data<-function(d){
  # Resample the data for training
  resampled_d = sample_n(d, nrow(d),  replace=T)
  # Keep original data as test and change all column names to test
  test_d = d
  colnames(test_d) = str_c('test', colnames(test_d), sep = '_')
  
  # Concatenate the datasets column wise.  Each has n observations, so we can use tidybayes::compose_data
  # To make a list of data to pass to stan
  stan_data = compose_data(bind_cols(resampled_d, test_d))
  
  return(stan_data)
}


extract_prediction<-function(fit, var_name){
  
  fit$draws(var_name) %>% 
    as_draws_matrix() %>% 
    apply(., 2, mean)
  
}


concentration = function(ts,C0, D, f, V, ka, ke){
  C0 + D*f*ka/(V*(ke-ka)) * (exp(-ka*ts) - exp(-ke*ts))
}



design_matrices = function(v.formula, ke.formula, f.formula, a.formula, data){
  stan_data = compose_data(data)
  stan_data$X_v = model.matrix(v.formula, data=d)
  stan_data$p_v = ncol(stan_data$X_v)
  
  stan_data$X_ke = model.matrix(ke.formula, data=d)
  stan_data$p_ke = ncol(stan_data$X_ke)
  
  stan_data$X_F = model.matrix(f.formula, data=d)
  stan_data$p_F = ncol(stan_data$X_F)
  
  stan_data$X_a = model.matrix(a.formula, data=d)
  stan_data$p_a = ncol(stan_data$X_a)
  
  return(stan_data)
}