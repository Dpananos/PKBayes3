library(tidyverse)
library(readxl)
library(janitor)
library(duckdb)
library(DBI)

con = dbConnect(duckdb(), dbdir='data/database/apixaban_data.duckdb')

# Write the raw data to the database for posterity
raw_data = read_xlsx('data/ute_full_raw_data.xlsx') %>% 
                clean_names()




# Write a cleaned version of the data, with some renamed variables
clean_1 = raw_data %>% 
  rename(age = age_at_enrollment,
         subjectids = new_subject_id,
         dose_mg_twice_daily = apixaban_dose_mg_twice_daily,
         yobs_ng_ml = apixaban_concentration_ng_ml,
         carbamazepine_mg_day = carbamazepine) %>% 
  select(
    subjectids,
    age, 
    sex, 
    weight_kg, 
    dose_mg_twice_daily,
    hrs_post_dose,
    yobs_ng_ml,
    creatinine_micromol_l,
    indication_for_doac_use,
    contains('hx'),
    contains('mg_day')
  ) %>% 
  mutate_at(vars(contains('hx')), list(~as.numeric(.=='yes')))

colnames(clean_1) = str_replace(colnames(clean_1), 'hx', 'history')


# Write a scaled version of the data for use in some modelling pursuits

clean_2 = clean_1 %>% 
  mutate(
    sex = if_else(sex=='male', 1, 0)
  ) %>% 
  rename(is_male=sex) %>% 
  mutate_at(vars(age, weight_kg, creatinine_micromol_l, contains('mg_day')), list(~as.vector(scale(.)))) %>% 
  select(-history_heart_failure, -ketoconazole_mg_day, -fluconazole_mg_day)
  



dbWriteTable(con, 'raw_data', raw_data, overwrite=F)
dbWriteTable(con, 'cleaned_data',clean_1, overwrite=F)
dbWriteTable(con, 'scaled_cleaned_data',clean_2, overwrite=F)

dbDisconnect(con)
