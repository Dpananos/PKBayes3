from matplotlib.pyplot import show
import pandas as pd
import numpy as np
from patients import Patient
import cmdstanpy
import json

def main():
    num_dense_patients = 40
    num_sparse_patients = 400

    # Create patients densley sampled.
    # See Patients.covars.  Sampled every half hour from 0.5 hours after ingestion to 12 hours after ingestion
    # Need to 'explode' by time because 'time' column is an array inside a column.
    dense_patients = [Patient(subjectid=i, sampling='dense').covars for i in range(num_dense_patients)]
    dense_df = pd.DataFrame(dense_patients).explode('time')

    # Create patients Sparsley sampled.
    # See Patients.covars.  Sampled sometime between 108 and 120 hours after ingestion.
    # This means hours post dose is time - 108 if doses are taken once every 12 hours with perfect adherence.
    sparse_patients = [Patient(subjectid=i+num_dense_patients, sampling='sparse').covars for i in range(num_sparse_patients)]
    sparse_df = pd.DataFrame(sparse_patients)

    # Combine the two simulated sets and recode their identifiers for passing to stan.
    df = pd.concat((dense_df, sparse_df))
    df['subjectids'] = pd.Categorical(df.subjectids).codes + 1

    # Create a dictionary to pass to the stan model
    stan_simulation_data = {
        'n': len(df),
        'n_subjectids': df.subjectids.nunique(),
        'subjectids': df.subjectids.tolist(),
        'time': df.time.tolist(),
        'age': df.age.tolist(),
        'weight': df.weight.tolist(),
        'creatinine': df.creatinine.tolist(),
        'is_male': df.is_male.tolist(),
        'dose': df.dose.tolist()
    }

    print(stan_simulation_data)



if __name__ == '__main__':

    with open('model_data.json') as file:
        model_data = json.loads(file.read())

        stan_file = '../models/008_combined_data_model.stan'
        model = cmdstanpy.CmdStanModel(stan_file=stan_file)
        fit = model.sample(model_data, show_progress=True, parallel_chains=4)
        fit.draws_pd().to_csv("draws.csv", index = False)