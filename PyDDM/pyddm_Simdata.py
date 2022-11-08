from ddm import Sample
from ddm.plot import model_gui
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm
from math import sin, radians, sqrt


# load in data
# for now, this is simulated data generated from Matlab code
# eventually, we can use Model below to simulate data too
df = pd.read_csv('/Users/stevenjerjian/Desktop/simKiani09_20200506.csv')

df["Trial"] = np.arange(len(df.index)) + 1 
# not sure why Unnamed is there, but ok..let's drop it
# and set the dataframe index as trial number

df = df.drop(['Unnamed: 7'], axis=1).set_index(['Trial'])

# create "correct" column
df["correct"] = ((df["choice"] == 2) & (df["heading"] > 0)) | ((df["choice"] == 1) & (df["heading"] < 0))
df["correct"] = df["correct"].astype(int)

# replace eps with zero (after assigning correct)
df["heading"].where(np.abs(df["heading"]) > 1e-3, 0, inplace=True)

# pyddm wants RT in seconds
df["RT"] /= 1000

# following the pyDDM example, use RT and correct as dep variables...should we be using choice instead?
sample = Sample.from_pandas_dataframe(df, rt_column_name="RT", correct_column_name="correct")

duration = 2

# define stimulus physical profiles
# Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
x = np.linspace(0, duration-0.001, duration*1000)
vel = norm.pdf(x, duration/2, 0.210)

vel = 0.37*(vel/max(vel))
acc = np.gradient(vel)*1000  # multiply by 1000 to get from m/s/ms to m/s/s

show_stim = False
if show_stim is True:
    fig = plt.figure(figsize=(10,10))
    plt.plot(x, vel)
    plt.plot(x, acc)
    plt.show()

vel = np.abs(vel/np.max(vel))
acc = np.abs(acc/np.max(acc))

# ~~~~~~~~~~~~~~~~~~~~~~#

# subclass drift to create a drift rate that varies according to dots3DMP conditions and stimulus physical profiles
import ddm.models

class DriftRate(ddm.models.Drift):
    name = '3DMP integrated drift rate'
    required_parameters = ['kves', 'kvis', 'acc', 'vel']
    required_conditions = ['heading', 'coherence', 'delta', 'modality']

    def get_drift(self, t, conditions, **kwargs):

        if conditions['modality'] == 1:
            return self.acc[int(t*1000)] * float(self.kves) * sin(radians(conditions['heading']))

        elif conditions['modality'] == 2:
            return self.vel[int(t*1000)] * float(self.kvis) * conditions['coherence'] * sin(radians(conditions['heading']))

        else:
            mu_ves = self.acc[int(t*1000)] * float(self.kves) \
                * sin(radians((conditions['heading'] - conditions['delta'] / 2)))

            mu_vis = self.vel[int(t*1000)] * float(self.kvis) * conditions['coherence'] \
                * sin(radians((conditions['heading'] + conditions['delta'] / 2)))

            wves = sqrt(float(self.kves)**2 / (float(self.kvis)**2 + float(self.kves)**2))
            wvis = sqrt(float(self.kvis)**2 / (float(self.kvis)**2 + float(self.kves)**2))

            return wves * mu_ves + wvis * mu_vis


from ddm import Model, Fittable
from ddm.functions import fit_adjust_model
from ddm.models import NoiseConstant, BoundConstant, OverlayNonDecision

model = Model(name="3DMP Kiani09 Sim",
              drift=DriftRate(kves=Fittable(minval=0, maxval=3), kvis=Fittable(minval=3, maxval=6), acc=acc, vel=vel),
              bound=BoundConstant(B=Fittable(minval=50, maxval=90)),
              noise=NoiseConstant(noise=1),
              overlay=OverlayNonDecision(nondectime=Fittable(minval=.1, maxval=.4)),
              dx=.01, dt=.001, T_dur=2.0)

fit_model = fit_adjust_model(sample=sample, model=model)

#model_gui(model=model, sample=sample)
