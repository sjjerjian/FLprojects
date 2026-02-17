# %% # DDM Testing Script


import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from behavior.Accumulator import Accumulator

from behavior.selfmotionddm import SelfMotionDDM, get_stim_urgs
from behavior.utils import dots3DMP_create_trial_list
from behavior.descriptive import plot_behavior_hdg, behavior_means, replicate_ves

def tic():
    return time.perf_counter()
def toc(tstart):
    return time.perf_counter() - tstart

# %% Initialize grid vectors for diffusion particle

# NOTE: 
# going too low on time_vec (like 0.01) can start to introduce some weirdness
# other values also have sweet spots e.g. bound ~1 and grid_vec -3 --> 0.
# e.g. try changing the bound to 3...weird stuff happens.

grid_vec = np.arange(-3, 0, 0.01) 
time_vec = np.arange(0, 2, 0.025)   

acc, vel = get_stim_urgs(time_vec)

# %% ==== demonstrate use of single Accumulator object ====

drifts = [0.25, 0.5, 1., 2.]
accum = Accumulator(
    grid_vec=grid_vec,
    tvec=time_vec,
    bound=1.5,
)
# this "applies" the drifts in [drifts] by creating [tvecx2] anti-correlated accumulators
accum.apply_drifts(drifts, labels=drifts) 

# this computes the cdf (for choice/RT) and pdf/log_odds (for wager) using method of images
tstart = tic()
accum.compute_distrs(return_pdf=True, use_vectorized=True) 
time_taken = toc(tstart)
print(f"Accumulator run took {time_taken:.3f} seconds")
accum.log_posterior_odds()
accum.plot();

# %% ==== Plot simulation of actual decision variable ====

# dotted lines show underlying accumulator drifts
# solid lines show simulated anti-correlated accumulators drawn from sampling
# black line shows the bound
dec_var, dv_fig = accum.dv(d_ind=1, show=True)

# %% ==== Initialize DDM object ====

init_params = {
    'kmult': [1, 1.5],        # ves, vis sensitivites
    'bound': [1, 1, 1],         # ves, vis, comb bounds
    'non_dec_time': [0.3],      # non-decision time (secs)
    'wager_thr': [0.5, 0.5, 0.5],     # log odds threshold for high bets
    'wager_alpha': [0.05],      # base rate of high bets
}

# generate list of unique conditions
# default delta=0 and nreps=1
X = dots3DMP_create_trial_list(
    hdgs=[-12, -6, -3, 0, 3, 6, 12],
    mods=[1, 2, 3],
    cohs=[0.3, 0.7],
)

# initialize DDM object (just for inference/prediction)
ddm = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,  # scale ves/vis according to acc/vel signals?
    return_wager=True   # whether to compute pdfs and log odds maps, and return wagers
    )

# %% ==== Generate model predictions ====

_, preds_model = ddm.predict(X, n_samples=5000, cache_accumulators=True)

# %% ==== Plot wager accumulator ====
ves_accum = ddm.accumulators_[('wager', 1.0)]
ves_accum.plot()


# %% ==== Simulate and visualize model predictions ====

# Note the errorbars here aren't "real" because we haven't really simulated trials,
# we've just drawn the probabilities for n_sample predictions for each of the unique conditions
# in general, we would use these preds_model as the data_fit points rather than the observations

preds_full = pd.concat((X, preds_model), axis=1)
preds_full = replicate_ves(preds_full) 
preds_full.sort_values(by=['modality', 'coherence', 'heading'])

df_means = behavior_means(preds_full, by_conds=['modality', 'coherence', 'heading'], long_format=True)

mod_map = {1: "ves", 2: "vis", 3: "comb"}
df_means['modality'] = df_means['modality'].map(mod_map)

plot_behavior_hdg(
    df_means,
    col='coherence',
    hue='modality',
    palette=['k', 'r', 'b'],
    hue_order=['ves', 'vis', 'comb'],
    )
# %%
