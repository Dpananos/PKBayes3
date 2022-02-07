import numpy as np
import pandas as pd
import duckdb
import pymc as pm
from sklearn.preprocessing import StandardScaler

con = duckdb.connect('data/database/apixaban_data.duckdb')

rdf = pd.read_sql('select * from rommel_cleaned_data', con=con)
udf = pd.read_sql('select * from ute_cleaned_data', con=con)

con.close()


# Get rommel's patient data

rommel_patients = rdf.drop_duplicates(['subjectids']).loc[:, ['age','weight_kg','creatinine_micromol_l']].values
ute_patients = udf.loc[:, ['age','weight_kg','creatinine_micromol_l']].values
covars = np.r_[rommel_patients, ute_patients]
scaler = StandardScaler().fit(covars)

rdf.loc[:, ['age','weight_kg','creatinine_micromol_l']] = scaler.transform(rdf.loc[:, ['age','weight_kg','creatinine_micromol_l']].values)
udf.loc[:, ['age','weight_kg','creatinine_micromol_l']] = scaler.transform(udf.loc[:, ['age','weight_kg','creatinine_micromol_l']].values)

# For modelling
r_ids = pd.Categorical(rdf.subjectids).codes
r_n_subjectids = np.unique(r_ids).size
r_age = rdf.age.values
r_is_male = (rdf.sex=='male').astype(int)
r_weight = rdf.weight_kg.values
r_creat = rdf.creatinine_micromol_l.values
r_time = rdf.hrs_post_dose.values
r_yobs = rdf.yobs_ng_ml.values*1000
r_D = rdf.dose_mg_twice_daily.values
X = np.c_[r_age, r_weight, r_creat, r_is_male]


u_ids = pd.Categorical(udf.subjectids).codes
u_age = udf.age.values
u_is_male = (udf.sex=='male').astype(int)
u_weight = udf.weight_kg.values
u_creat = udf.creatinine_micromol_l.values
u_time = udf.hrs_post_dose.values
u_yobs = udf.yobs_ng_ml.values*1000
u_D = udf.dose_mg_twice_daily.values
u_amio = udf.amiodarone_mg_day.values


def concentration(time, dose, f, cl, ke, ka, c0=0.0):

    kernel = pm.math.exp(-ka*time) - pm.math.exp(-ke*time)
    factor = c0 + dose*f*ke*ka / (cl * (ke - ka)) 
    return kernel*factor
    

nrows, ncols = X.shape
with pm.Model() as model:

    # Pharmacokientic parameters

    mu_cl = pm.Normal('mu_cl', 1.19, 0.1)
    mu_t = pm.Normal('mu_t', 1.19, 0.1)
    mu_a = pm.Normal('mu_a', -0.25, 0.5)
    mu_F = pm.Normal('mu_F', 0,0.025)

    beta_cl = pm.Normal('beta_cl', 0, 0.25, shape = ncols)
    beta_t = pm.Normal('beta_t', 0, 0.25, shape = ncols)
    beta_a = pm.Normal('beta_a', 0, 0.25, shape = ncols)

    s_cl =  pm.Gamma('s_cl', 15, 100)
    s_t =  pm.Gamma('s_t', 5, 100)
    s_a =  pm.Gamma('s_a', 10, 100)

    z_cl = pm.Normal('z_cl', 0, 1, shape=r_n_subjectids)
    z_t = pm.Normal('z_t', 0, 1, shape=r_n_subjectids)
    z_a = pm.Normal('z_a', 0, 1, shape=r_n_subjectids)

    cl_linpred = mu_cl + pm.math.dot(X, beta_cl) + s_cl*z_cl[r_ids]
    t_linpred = mu_t + pm.math.dot(X, beta_t) + s_t*z_t[r_ids]
    a_linpred = mu_a + pm.math.dot(X, beta_a) + s_a*z_a[r_ids]

    cl = pm.Deterministic('cl', pm.math.exp(cl_linpred))
    tmax = pm.Deterministic('tmax', pm.math.exp(t_linpred))
    a = pm.Deterministic('a', pm.math.invlogit(a_linpred))
    F = pm.Deterministic('F', pm.math.invlogit(mu_F))
    ka = pm.Deterministic('ka', pm.math.log(a) / (tmax * (a-1)))
    ke = pm.Deterministic('ke', pm.math.log(a) * a/ (tmax * (a-1)))

    # Time delay
    phi = pm.Beta('phi', 20, 20)
    kappa = pm.Beta('k', 20, 20)
    delta = pm.Beta('delta', phi/kappa, (1-phi)/kappa, shape=r_n_subjectids)
    delayed_time = r_time - 0.5*delta[r_ids]

    conc  = concentration(delayed_time, r_D, F, cl, ke, ka)
    latent_conc = pm.Deterministic('latent_conc', conc)

    sigma = pm.Lognormal('sigma', -2.3, 0.2)
    Y = pm.Lognormal('Y', pm.math.log(conc), sigma=sigma,  observed=r_yobs)


    trace = pm.sample(chains=4)