library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(rms)
library(posterior)
library(DBI)
library(duckdb())
source('R/utils.R')

theme_set(theme_light())

con = dbConnect(duckdb(), 'data/database/apixaban_data.duckdb')
d = tbl(con, "cleaned_Data") %>% 
    collect()
dbDisconnect(con)

dd = datadist(d)
options(datadist='dd')
markus_model = ols(log(yobs_ng_ml) ~ dose_mg_twice_daily + hrs_post_dose + sex + age + weight_kg + creatinine_micromol_l + amiodarone_mg_day + diltiazem_mg_day, data = d, x=T, y=T)
extended_model = ols(log(yobs_ng_ml) ~ dose_mg_twice_daily + hrs_post_dose*(sex + age + weight_kg + creatinine_micromol_l) + amiodarone_mg_day + diltiazem_mg_day, data = d, x=T, y=T)
model_005_tmax_priors = cmdstan_model('models/model_005_tmax_priors.stan')

# Validate each model using optimism corrected bootstrap
validate(markus_model, B = 250)
validate(extended_model, B = 250)

bootstrapped_optimism = replicate(20, {
  
  stan_data = prepare_bayesian_data(d)
  fit = model_005_tmax_priors$sample(stan_data, chains=4, parallel_chains=4)
  
  # Test performance
  ypred = extract_prediction(fit, 'test_yppc')
  ytest = stan_data$test_yobs_ng_ml
  test_performance = MLmetrics::MSE(log(ypred), log(ytest))
  
  # Bootstrap performance
  ypred = extract_prediction(fit, 'yppc')
  ytrue = stan_data$yobs_ng_ml
  bootstrap_performance = MLmetrics:: MSE(log(ypred), log(ytrue))
  
  optimism = bootstrap_performance - test_performance
  
})

print(apparent_performance - mean(bootstrapped_optimism))

