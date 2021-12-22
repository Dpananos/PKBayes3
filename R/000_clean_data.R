library(tidyverse)
library(readxl)
library(janitor)
library(duckdb)
library(DBI)

raw_data = read_xlsx('data/ute_full_raw_data.xlsx') %>% 
                clean_names()



con = dbConnect(duckdb(), dbdir='data/database/apixaban_data.duckdb')

dbWriteTable(con, 'raw_data', raw_data, overwrite=T)

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
  )

colnames(clean_1) = str_replace(colnames(clean_1), 'hx', 'history')

dbWriteTable(con, 'cleaned_data',clean_1, overwrite=T)

dbDisconnect(con)
