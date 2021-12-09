library(tidyverse)
library(readxl)
library(janitor)


raw_data = read_xlsx('data/ute_full_raw_data.xlsx') %>% 
           clean_names()



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

clean_1 %>% 
  mutate_at(vars(contains('history')), ~as.numeric(.=='Yes')) %>% 
  write_csv('data/cleaned_data.csv')
