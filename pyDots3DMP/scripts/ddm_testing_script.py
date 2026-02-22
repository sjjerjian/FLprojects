# %% ================================================
# DDM DEMO SCRIPT
# ===================================================

import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

# set a save location
save_dir = "param_recov"
Path.mkdir(Path(save_dir), parents=True, exist_ok=True)

# Set up logger to console and file - this will show us the logging info for our DDM, and the bads fitting routine
logger = logging.getLogger()
logger.setLevel(logging.INFO)
fmt = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    log_path = Path(__file__).resolve().parent / "ddm_testing.log"
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

# custom imports
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
drifts = [0, 0.25, 0.5, 1., 2.]
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
# kmult/bound will be used to determine drift rates/bounds for underlying Accumulator objects (1 per condition)
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
ddm_obj = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,  # scale ves/vis according to acc/vel signals?
    return_wager=True   # whether to compute pdfs and log odds maps, and return wagers
    )

# save to disk, reload from disk, compare model parameters
# ddm.save("model.json")
# ddm2 = SelfMotionDDM.load("model.json")
# print(SelfMotionDDM.params_table(ddm, ddm2))

# here we generate 100 trials of each condition, but set n_samples to 1 to get model predictions
# of choice and wager as binary outcomes for each trial
X = dots3DMP_create_trial_list(
    hdgs=[-12, -6, -3, 0, 3, 6, 12],
    mods=[1, 2, 3],
    cohs=[0.3, 0.7],
    nreps=100,
)

# predict returns model_probs, sampled_predictions
_, preds_model = ddm_obj.predict(X, n_samples=1, cache_accumulators=True)

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

# plot simulated data. currently no fit curve, so will just draw lines between each
# can fit with a gaussian eventually (using behavior.utils.gauss_fit_hdg_group)
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

# start from a few different points to what we used to generate the data, see if the model can recover
init_params2 = {
    'kmult': [1.5, 1.0],            # ves, vis sensitivites
    'bound': [0.3, 0.6, 1.2],       # ves, vis, comb bounds
    'non_dec_time': [0.6],          # non-decision time (secs)
    'wager_thr': [0.5],   # log odds threshold for high bets
    'wager_alpha': [0.05],          # base rate of high bets
}

# initialize new DDM object
ddm_fit = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params2, 
    stim_scaling=True,  
    return_wager=True,
    save_dir=save_dir
    )

# set some options for the BADS routine (or scipy.minimize)
fit_options = {
        "random_seed": 42,
        "max_fun_evals": 100,
        "display": "full"
    }

# run the fit, with some fixed params
ddm_fit.fit(
    X,
    preds_model,   # simulated choice, PDW, RT from initial setup
    fixed_params=["wager_alpha"],  
    fit_method='bads',
    fit_options=fit_options
)

# print comparison table of params
# TODO flip rows and columns here, or don't bother with in-built method...
print(SelfMotionDDM.params_table(ddm_obj, ddm_fit))


# %% ================================================
# Plot fitted curves on top of original simulated data
# # ===================================================

# create a nicely spaced set of headings to run predictions over using fitted parameters
X_pred = dots3DMP_create_trial_list(
    hdgs=np.linspace(-12, 12, 100),
    mods=[1, 2, 3],
    cohs=[0.3, 0.7],
    nreps=1,        # only need 1
    shuff=False,
)

# use the wager maps from fitting to the short list of headings, i.e. DON'T recompute them with the new heading set!!
preds_, preds_samples = ddm_fit.predict(X_pred, n_samples=1, use_cached_wager_maps=True, rt_sampling_method="mean")

# a couple of kluges here to get the nice model predictions dataframe
# RT predictions come from the samples df, preds_ originally has the RT likelihoods needed for fitting
# for choice and PDW we can use preds_ columns as is, because they are probabilities of the binary outcome

preds_['RT'] = preds_samples['RT']
preds_full = pd.concat((X_pred, preds_), axis=1)
preds_full = replicate_ves(preds_full) 
mod_map = {1: "ves", 2: "vis", 3: "comb"}
preds_full['modality'] = preds_full['modality'].map(mod_map)

g = plot_behavior_hdg(
    df_means,
    data_fit=preds_full,
    col='coherence',
    hue='modality',
    palette=['k', 'r', 'b'],
    hue_order=['ves', 'vis', 'comb'],
    )
g.figure.savefig(save_dir / "behavior_fit.png", dpi=150, bbox_inches="tight")
# %%
