import cmdstanpy
import json
import os



with open('model_data.json') as file:
    model_data = json.loads(file.read())
    model = cmdstanpy.CmdStanModel(stan_file='../models/008_combined_data_model.stan')
    fit = model.sample(model_data, chains=4, parallel_chains=4, show_progress=True)