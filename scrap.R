library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(posterior)
library(duckdb)
library(DBI)
source('R/utils.R')

theme_set(theme_light())


# Make a connection to the database containing the data for the study.
# This command allows R to talk to the database
con = dbConnect(duckdb(), 'data/database/apixaban_data.duckdb')

# I've saved the data into a database rather than csvs.
# csvs are fine, but there is a lot of clutter that goes on when you create many derived files.
# Extract the data using tools from dbplyr
# The result is a tibble containing the data.
d = tbl(con, "cleaned_Data") %>% 
    collect() %>% 
    mutate(i = seq_along(subjectids)) %>% 
    select(-history_heart_failure)

# Close the connection.  This is good code hygiene.
# If you want to edit the d variable (e.g. ignoring a column, like I have done), you need to first reconnect (line 14)
dbDisconnect(con, shutdown=TRUE)

# Stan requires the data to be passed to our model in a specific way.
# The desired format is a list of vectors.
# tidybayes::compose_data is a handy function for taking tidy data
# and turning it into a Stan friendly form.
# There are some sharp corners, I encourage you to read the docs.
stan_data = compose_data(d)

# Load the stan model and sample from it.
model = cmdstan_model('models/scrap_models.stan')
fit = model$sample(stan_data, chains=4, parallel_chains=4)

#----

# Predicted versus truth
fit %>% 
  spread_draws(yppc[i], ndraws = 250) %>% 
  mean_qi %>% 
  left_join(d) %>% 
  ggplot(aes(yppc,yobs_ng_ml))+
  geom_abline()+
  geom_point(aes(color = factor(dose_mg_twice_daily)))+
  geom_smooth()


curves = fit %>% 
  spread_draws(
    C0[i],
    f[i],
    V[i],
    ke[i], 
    ka[i]
  ) %>% 
  left_join(select(d, i, dose_mg_twice_daily)) %>% 
  crossing(ts = seq(0, 12)) %>% 
  mutate(y = concentration(ts,C0, dose_mg_twice_daily, f, V, ka, ke)) %>% 
  group_by(i, ts) %>% 
  mean_qi(y)

curves %>% 
  left_join(d) %>% 
  ggplot(aes(ts, 1000*y, group = subjectids))+
  geom_line(aes(color = factor(sex)),  alpha = 0.5)+
  facet_grid(sex~dose_mg_twice_daily)+
  geom_point(data = d, aes(hrs_post_dose, yobs_ng_ml, fill = factor(sex)), shape=21, color='white')+
  labs(x='Hrs Post Dose', y='Concentration (ng/ml)')



d$ke = extract_prediction(fit, 'ke')


d %>% 
  ggplot(aes(ke, sex))+
  geom_point()
