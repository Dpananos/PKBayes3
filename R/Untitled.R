library(tidyverse)
library(duckdb)
library(DBI)
library(cmdstanr)
library(tidybayes)
library(posterior)
source('R/utils.R')
# The point of this exercise is to use Romme's very high quality data to 
# calibrate our model in some sense.
# The PK should be shared in some way between the two populations.
# Obviously, there will be some variability, but our hope is that there won't be too much and we can pool what we learn
# from Rommel's data and apply it to Ute's data.
# Let's set up some data to pass to a Stan model.


# Query data from the database
con = dbConnect(duckdb(),'data/database/apixaban_data.duckdb')

ute_data = tbl(con, 'ute_cleaned_data') %>% 
           collect()

rommel_data = tbl(con, 'rommel_cleaned_data') %>% 
              collect()

dbDisconnect(con)


# Prepare data to pass to stan
# Rommel first
distinct_r = distinct(rommel_data, subjectids, .keep_all = T) %>% 
             mutate(i = seq_along(subjectids))


r_n = nrow(rommel_data)
r_subjectids = rommel_data$subjectids
r_n_subjectids = n_distinct(r_subjectids)
r_time = rommel_data$hrs_post_dose
r_yobs = rommel_data$yobs_ng_ml

r_age  = distinct_r$age
r_weight  = distinct_r$weight_kg
r_creatinine  = distinct_r$creatinine_micromol_l
r_is_male = if_else(distinct_r$sex=='male', 1, 0)
r_D = distinct_r$dose_mg_twice_daily

r_data = list(
  r_n = r_n,
  r_n_subjectids = r_n_subjectids,
  r_subjectids = r_subjectids,
  r_time = r_time,
  r_yobs = r_yobs,
  r_age = r_age,
  r_weight = r_weight,
  r_creatinine = r_creatinine,
  r_is_male = r_is_male,
  r_D = r_D
)

# Now Ute


u_time = ute_data$hrs_post_dose
u_yobs = ute_data$yobs_ng_ml
u_age  = ute_data$age
u_weight  = ute_data$weight_kg
u_creatinine  = ute_data$creatinine_micromol_l
u_is_male = if_else(ute_data$sex=='male', 1, 0)
u_n = nrow(ute_data)
u_D = ute_data$dose_mg_twice_daily
u_amiodarone = ute_data$amiodarone_mg_day
u_diltiazem = ute_data$diltiazem_mg_day

u_data = list(
  u_n = u_n,
  u_time = u_time,
  u_yobs = u_yobs,
  u_age = u_age,
  u_weight = u_weight,
  u_creatinine = u_creatinine,
  u_is_male = u_is_male,
  u_D = u_D,
  u_amiodarone=u_amiodarone,
  u_diltiazem=u_diltiazem
)

model_data = c(r_data, u_data)

# Load model  
if(T){
  model = cmdstan_model('models/008_combined_data_model.stan')
  fit = model$sample(model_data, chains=4, seed=0, parallel_chains=4)
  
  fit$save_object(file='model_008.RDS')
  
}



fit %>% 
  spread_draws(z_F[i]) %>% 
  mean_qi %>% 
  bind_cols(distinct(rommel_data, subjectids, .keep_all = T)) %>% 
  ggplot(aes(weight_kg, z_F, ymin=.lower, ymax=.upper))+
  geom_pointrange()

ypred = extract_prediction(fit, 'u_C')*1000
ytrue = u_data$u_yobs
plot(log(ypred), log(ytrue))
abline(0, 1)
MLmetrics::RMSE(log(ypred), log(ytrue))


fit %>% 
  spread_draws(r_yppc[i]) %>% 
  mean_qi %>% 
  bind_cols(rommel_data) %>% 
  ggplot(aes(hrs_post_dose, 1000*r_yppc, ymin = .lower*1000, ymax = .upper*1000))+
  geom_line()+
  geom_point(data=rommel_data, aes(hrs_post_dose, yobs_ng_ml), inherit.aes = F)+
  geom_ribbon(alpha=0.5)+
  facet_wrap(~subjectids)


