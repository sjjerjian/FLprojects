# DDM Testing Script

# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from behavior.Accumulator import Accumulator

from behavior.selfmotionddm import SelfMotionDDM
from behavior.utils import dots3DMP_create_trial_list
from behavior.descriptive import plot_behavior_hdg, gauss_fit_hdg_group


# %% Initialize grid vectors for diffusion particle
grid_vec = np.arange(-3, 0, 0.05)  # -3 to 0, steps of 0.01
time_vec = np.arange(0, 2, 0.01)   # 0 to 2s, steps of 50ms

# %% demonstrate use of single Accumulator object

drifts = [0, 0.25, 0.5, 1., 2.]
accum = Accumulator(
    grid_vec=grid_vec,
    tvec=time_vec,
    bound=1,
)
accum.apply_drifts(drifts, labels=drifts) # this "applies" the drifts in [drifts] by creating [tvecx2] anti-correlated accumulators
accum.compute_distrs(return_pdf=True, use_vectorized=True) # this computes the cdf (for choice/RT) and pdf (for wager) using method of images
#accum.log_posterior_odds()

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
    return_wager=True       # whether to compute pdfs and log odds maps, and return wagers
    )

preds = ddm.simulate(X, n_samples=5, sample_dvs=True, seed=1)

# %%
