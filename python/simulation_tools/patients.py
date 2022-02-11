import numpy as np
import random
from scipy.stats import norm, binom, uniform


class BasePatient:

    def __init__(self):

        self._age_mean = 50
        self._age_std = 11

        self._weight_mean = 88
        self._weight_std = 24

        self._creatine_mean = 68
        self._creatinine_sd = 12

        self._doses = [2.5, 5.0]


    def draw_covar(self, mu, sd, round_to=2):

        covar = norm(loc=mu, scale=sd).rvs()
        
        return np.round(covar, round_to)


class Patient(BasePatient):

    def __init__(self, subjectid, sampling=None):

        super().__init__()
        # Draw the z score for easier passing to stan
        # We can get the raw result by adding mean and sd.
        self.subjectid = f'subejct_{subjectid:03}'
        self.age = self.draw_covar(0, 1, round_to=1)
        self.weight = self.draw_covar(0, 1, round_to=2)
        self.creatinine = self.draw_covar(0, 1, round_to=2)
        self.is_male = random.choice([0.0, 1.0])
        self.dose = random.choice([2.5, 5.0])
        self.sampling = sampling

        # Assume much like our real data that our controlled study has no one taking amiodarone,
        # leaving us to estimate the effect from the sparse data.
        if self.sampling == 'sparse':
            self.amio = random.choice([0, 0.25, 0.5, 0.75, 1])
        else:
            self.amio = 0

    @property
    def covars(self):

        d = {
            'subjectids': self.subjectid,
            'age': self.age,
            'weight': self.weight,
            'creatinine': self.creatinine,
            'is_male': self.is_male,
            'dose': self.dose,
            'amio': self.amio,
            'sampling': self.sampling
        }

        if self.sampling == 'dense':
            t = np.arange(start=0.5, stop=12.5, step=0.5)
        elif self.sampling == 'sparse':
            t = np.round(uniform(108, 12).rvs(), 2)
        else:
            raise ValueError('sampling must be "dense" or "sparse".')
        
        d['time'] = t

        return d
