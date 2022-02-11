import numpy as np
import pandas as pd


def create_simulation_stan_data(df, amio_effect):

    stan_simulation_data = {
        'n': len(df),
        'n_subjectids': df.subjectids.nunique(),
        'subjectids': df.subjectids.tolist(),
        'time': df.time.tolist(),
        'age': df.age.tolist(),
        'weight': df.weight.tolist(),
        'creatinine': df.creatinine.tolist(),
        'is_male': df.is_male.tolist(),
        'dose': df.dose.tolist(),
        'amio': df.amio.tolist(),
        'amio_effect': amio_effect
    }

    return stan_simulation_data



def prepare_data(df, prefix):


    patients = df.drop_duplicates(['subjectids'])
    
    # Need to re-encoode the ids since I am passing dense/sparse data seperately.
    # Ids must be unique to each set.
    df['subjectids'] = pd.Categorical(df.subjectids).codes + 1

    simulation_covars = {
        f'{prefix}_age': patients.age.tolist(),
        f'{prefix}_weight': patients.weight.tolist(),
        f'{prefix}_creatinine': patients.creatinine.tolist(),
        f'{prefix}_is_male': patients.is_male.tolist(),
        f'{prefix}_dose': patients.dose.tolist(),
        f'{prefix}_amio': patients.amio.tolist()
    }

    simulation_data = {
        f'{prefix}_n': len(df),
        f'{prefix}_subjectids': df.subjectids.tolist(),
        f'{prefix}_n_subjectids': df.subjectids.nunique(),
        f'{prefix}_time': df.time.tolist(),
        f'{prefix}_yobs': df.yobs.tolist()
    }

    data = {**simulation_data, **simulation_covars}

    return data


def prepare_simulation_data(sparse_df, dense_df):

    dense_sim_data = prepare_data(dense_df, prefix='dense')
    sparse_sim_data = prepare_data(sparse_df, prefix='sparse')

    return {**dense_sim_data, **sparse_sim_data}



