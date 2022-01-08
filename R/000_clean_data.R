library(tidyverse)
library(readxl)
library(janitor)
library(duckdb)
library(DBI)

con = dbConnect(duckdb(), dbdir='data/database/apixaban_data.duckdb')

# Write the raw data to the database for posterity
raw_data = read_xlsx('data/ute_full_raw_data.xlsx') %>% 
                clean_names()

rommel_data = read_csv("data/rommel_data.csv")


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
  

# Clean rommel's data to match ute's

rommel_clean = rommel_data %>% 
  mutate(dose_mg_twice_daily = 2.5,
         Sex = str_to_lower(Sex),
         Subject = as.numeric(as.factor(Subject))) %>% 
  rename(
    hrs_post_dose = Time,
    subjectids = Subject,
    yobs_ng_ml = Concentration,
    age = Age,
    weight_kg = Weight,
    creatinine_micromol_l = Creatinine,
    sex = Sex
  ) %>% 
  select(
    subjectids,
    dose_mg_twice_daily,
    hrs_post_dose, 
    yobs_ng_ml,
    age,
    weight_kg, 
    creatinine_micromol_l,
    sex
  )


dbWriteTable(con, 'ute_raw_data', raw_data, overwrite=T)
dbWriteTable(con, 'ute_cleaned_data',clean_1, overwrite=T)
dbWriteTable(con, 'ute_scaled_cleaned_data',clean_2, overwrite=T)
dbWriteTable(con, 'rommel_cleaned_data',rommel_clean, overwrite=T)
dbWriteTable(con, 'rommel_raw',rommel_data, overwrite=T)

dbDisconnect(con)
