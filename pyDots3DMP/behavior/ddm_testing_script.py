# DDM Testing Script

# %%
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from behavior.Accumulator import Accumulator

from behavior.selfmotionddm import SelfMotionDDM, get_stim_urgs
from behavior.utils import dots3DMP_create_trial_list
from behavior.descriptive import plot_behavior_hdg, gauss_fit_hdg_group

def tic():
    return time.perf_counter()
def toc(tstart):
    return time.perf_counter() - tstart

# %% Initialize grid vectors for diffusion particle

# NOTE: going too low on time_vec (like 0.01) can start to introduce some weirdness
grid_vec = np.arange(-3, 0, 0.005) 
time_vec = np.arange(0, 2, 0.025)   

acc, vel = get_stim_urgs(time_vec)

# %% demonstrate use of single Accumulator object

drifts = [0, 0.25, 0.5, 1., 2.]
accum = Accumulator(
    grid_vec=grid_vec,
    tvec=time_vec,
    bound=1,
)
accum.apply_drifts(drifts, labels=drifts) # this "applies" the drifts in [drifts] by creating [tvecx2] anti-correlated accumulators

tstart = tic()
accum.compute_distrs(return_pdf=False, use_vectorized=True) # this computes the cdf (for choice/RT) and pdf (for wager) using method of images
time_taken = toc(tstart)
print(f"Accumulator run took {time_taken:.3f} seconds")
#accum.log_posterior_odds()
accum.plot();

# %%

dec_var, dv_fig = accum.dv(d_ind=1, show=True)

# %%

init_params = {
    'kmult': [1, 1],        # ves, vis sensitivites
    'bound': [1, 1, 1],         # ves, vis, comb bounds
    'non_dec_time': [0.3],      # non-decision time (secs)
    'wager_thr': [0.5, 0.5, 0.5],     # log odds threshold for high bets
    'wager_alpha': [0.05],      # base rate of high bets
}

# generate list of unique conditions
X = dots3DMP_create_trial_list(
    hdgs=[-12, -6, -3, 0, 3, 6, 12],
    mods=[1, 2, 3],
    cohs=[0.3, 0.7],
    deltas=[0],
    nreps=1
)

# initialize DDM object (just for inference/prediction)
ddm = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,      # scale ves/vis according to acc/vel signals? can also pass in a 2-length tuple of arrays
    return_wager=False       # whether to compute pdfs and log odds maps, and return wagers
    )

# %%

_, preds_model = ddm.predict(X, n_samples=1, cache_accumulators=True)

# %%
ves_accum = ddm.accumulators_[(1.0, 0.3, 0.0)]
ves_accum.plot()

# %%
# preds_simul = ddm.simulate(X, n_samples=5, sample_dvs=True, seed=1)
all_data = pd.concat((X, preds_model), axis=1)
all_data.sort_values(by=['modality', 'coherence', 'heading'])
# %%
