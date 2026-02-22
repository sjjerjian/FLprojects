# %% # DDM Demo Script

import time

import numpy as np
import pandas as pd

from ddm import Accumulator, SelfMotionDDM
from behavior.utils import dots3DMP_create_trial_list
from behavior.descriptive import plot_behavior_hdg, behavior_means, replicate_ves

# %% ================================================
# demonstrate use of Accumulator
# ===================================================

# NOTE: 
# going too low on time_vec (like 0.01) can start to introduce some weirdness
# other values also have sweet spots e.g. bound ~1 and grid_vec -3 --> 0.
grid_vec = np.arange(-3, 0, 0.01) 
time_vec = np.arange(0, 2, 0.025)  

# instantiate an Accumulator object
# (sets up the MOI quadrants and initial particle location)
accum = Accumulator(
    grid_vec=grid_vec,
    tvec=time_vec,
    bound=1.5,
)

# this "applies" the drifts in [drifts] by creating [tvecx2] anti-correlated accumulators
drifts = [0.25, 0.5, 1., 2.]
accum.apply_drifts(drifts, labels=drifts) 

# this computes the cdf (for choice/RT) and pdf/log_odds (for wager) using method of images
tstart = time.perf_counter()
accum.compute_distrs(return_pdf=True, use_vectorized=True) 
print(f"Accumulator run took {time.perf_counter() - tstart:.3f} seconds")
accum.log_posterior_odds()
# accum.plot();

# %% ================================================
# plot simulation of actual decision variable
# # ===================================================

# dotted lines show underlying accumulator drifts
# solid lines show simulated anti-correlated accumulators drawn from sampling
# black line shows the bound

dec_var, dv_fig = accum.dv(d_ind=1, show=True)

# %% ================================================
# Use the SelfMotionDDM setup for ves-vis cue comb task
# # ===================================================

# set initial parameters
# kmult/bound will be used to determine drift rates/bounds for underlying Accumulator objects
# other parameters are post-accumulator results, for yielding final behavioral variables

# for kmult, a list of length 2 implies that vis kmult will be scaled by coherence
# a list of length 3 gives low and high vis independent sensitivities
# for all parameters except kmult, a list of length 1 enforces the same value for all modalities,
# whereas a list of length 3 gives separate parameters to each modality

init_params = {
    'kmult': [1, 1.5],          # ves, vis sensitivites. if length 2, vis will be scaled by coh
    'bound': [1],               # ves, vis, comb bounds. 
    'non_dec_time': [0.3],      # non-decision time (secs)
    'wager_thr': [0.8],         # log odds threshold for high bets
    'wager_alpha': [0],      # base rate of high bets
}

# initialize DDM object
ddm = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,  # scale ves/vis according to acc/vel signals?
    return_wager=True   # whether to compute pdfs and log odds maps, and return wagers
    )

# here we generate 100 trials of each condition, but set n_samples to 1 to get model predictions
# of choice and wager as binary outcomes for each trial
X = dots3DMP_create_trial_list(
    hdgs=[-12, -6, -3, 0, 3, 6, 12],
    mods=[1, 2, 3],
    cohs=[0.3, 0.7],
    nreps=100,
)

# predict returns model_probs, sampled_predictions
_, preds_model = ddm.predict(X, n_samples=1, cache_accumulators=True)

# %% ================================================
# Visualize model 'predictions' i.e. simulated data
# # ===================================================

preds_full = pd.concat((X, preds_model), axis=1)
preds_full = replicate_ves(preds_full) # replicate ves for high coherence, for plotting convenience
preds_full.sort_values(by=['modality', 'coherence', 'heading'])

# plot_behavior_hdg expects means
df_means = behavior_means(preds_full, by_conds=['modality', 'coherence', 'heading'], long_format=True)

# map ordinal modality labels to strings
mod_map = {1: "ves", 2: "vis", 3: "comb"}
df_means['modality'] = df_means['modality'].map(mod_map)

# plot_behavior_hdg(
#     df_means,
#     col='coherence',
#     hue='modality',
#     palette=['k', 'r', 'b'],
#     hue_order=['ves', 'vis', 'comb'],
#     )

# %% ================================================
# Smoke test the model - start from new params and try and recover the original params
# # ===================================================

# for reference
# init_params = {
#     'kmult': [1, 1.5],          # ves, vis sensitivites. if length 2, vis will be scaled by coh
#     'bound': [1],               # ves, vis, comb bounds. 
#     'non_dec_time': [0.3],      # non-decision time (secs)
#     'wager_thr': [0.8],         # log odds threshold for high bets
#     'wager_alpha': [0.05],      # base rate of high bets
# }

# start from a few different points to before
init_params2 = {
    'kmult': [1.5, 1.0],            # ves, vis sensitivites
    'bound': [0.5, 0.5, 0.5],       # ves, vis, comb bounds
    'non_dec_time': [0.3],          # non-decision time (secs)
    'wager_thr': [0.5],   # log odds threshold for high bets
    'wager_alpha': [0],          # base rate of high bets
}

# initialize new DDM object
ddm_fit = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,  
    return_wager=True,
    )

fit_options = {
        "random_seed": 42,
        "max_fun_evals": 20,
        "display": "full"
    }

ddm_fit.fit(
    X,
    preds_model,   # simulated choice, PDW, RT from initial setup
    fixed_params=["non_dec_time", "wager_alpha"],   # fix these to reduce complexity
    fit_method='bads',
    fit_options=fit_options
)
print(ddm_fit.params_)

# %% ================================================
# Plot fitted curves on top of original simulated data
# # ===================================================

X_pred = dots3DMP_create_trial_list(
    hdgs=list(range(-12, 12)),
    mods=[1, 2, 3],
    cohs=[0.3, 0.7],
    nreps=1,
    shuff=False,
)

preds_['RT'] = preds_samples['RT']

preds_full = pd.concat((X_pred, preds_), axis=1)
preds_full = replicate_ves(preds_full) 
mod_map = {1: "ves", 2: "vis", 3: "comb"}
preds_full['modality'] = preds_full['modality'].map(mod_map)

plot_behavior_hdg(
    df_means,
    data_fit=preds_full,
    col='coherence',
    hue='modality',
    palette=['k', 'r', 'b'],
    hue_order=['ves', 'vis', 'comb'],
    )
# %%
