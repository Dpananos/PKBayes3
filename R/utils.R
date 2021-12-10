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
