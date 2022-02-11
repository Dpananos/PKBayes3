import pandas as pd
import numpy as np
import cmdstanpy
from itertools import product

from simulation_tools.patients import Patient
from simulation_tools.data_wrangling import create_simulation_stan_data, prepare_simulation_data

def main(num_dense_patients, num_sparse_patients,run, amio_effect):

    # Create patients densely sampled.
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
    stan_simulation_data = create_simulation_stan_data(df, amio_effect)

    # Load up simulation model and draw parameters
    simulation_model =cmdstanpy.CmdStanModel(stan_file='simulation_models/simulate_patients.stan')
    fit = simulation_model.sample(stan_simulation_data, chains=1, iter_warmup=0, iter_sampling=1, fixed_param=True )


    # Append simulation results to the data
    df['yobs'] = fit.stan_variable('observed_concentration').ravel()
    df['latent_y'] = fit.stan_variable('concentration').ravel()


    simulation_data = prepare_simulation_data(
        sparse_df=df.query('sampling=="sparse"'),
        dense_df = df.query('sampling=="dense"')
    )

    inference_model = cmdstanpy.CmdStanModel(stan_file='simulation_models/simulation_inference.stan')
    fit = inference_model.sample(simulation_data, parallel_chains=4, show_progress=True)

    draws = fit.draws_pd()

    estimated_amio_effect = draws[['b_F_amio']].describe(percentiles = [0.025, 0.05, 0.5, 0.95, 0.975])

    # Keep the index, it tells us what each statistic is
    estimated_amio_effect.to_csv(f'simulation_data/amio_effect_{amio_effect}_num_dense_patients_{num_dense_patients}_run_{run}.csv', index = True)






if __name__ == '__main__':

    amio_effects = [0, 0.125, 0.25, 0.5, 1.0, 1.5]
    nruns = np.arnage(1, 11)
    num_dense_patients = [10, 20, 30, 40, 50]
    p = product(num_dense_patients, nruns, amio_effects)

    for dense_patients, run, amio_effect in p:
        main(dense_patients, 10*dense_patients, run, amio_effect)
