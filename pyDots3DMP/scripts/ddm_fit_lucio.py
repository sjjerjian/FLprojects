"""Fit some real data"""

# TODO hand-fit as well as possible?
# add loop over multiple starting points / same starting point w different seeds
# add option to set bounds without params themselves

from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd

from ddm import SelfMotionDDM
from behavior.utils import dots3DMP_create_trial_list, data_cleanup
from utils import setup_loggers
import behavior.descriptive as behav

# %% ================================================ 
# set up save location and logger, and load data
# ===================================================
save_dir = Path(f"results/lucio/{datetime.now().strftime('%y%m%d_%H%M%S')}")
Path.mkdir(save_dir, parents=True, exist_ok=True)
setup_loggers(log_file_path=save_dir / "fit_lucio.log")

datafilepath = "/Users/stevenjerjian/FLprojects/lucio_20220512-20230606.csv"
data = data_cleanup(datafilepath)  # kluge to load and clean raw data, could be improved

# %% ================================================
# set up for fitting
# ===================================================

def process_predictions(X, y=None):
    """replicate ves for high coherence, for plotting convenience"""
    mod_map = {1: "ves", 2: "vis", 3: "comb"}
    data = pd.concat((X, y), axis=1) if y is not None else X.copy()
    data = behav.replicate_ves(data)
    data['modality'] = data['modality'].map(mod_map)
    return data

data_proc = process_predictions(data)

# plot_behavior_hdg function expects pre-computed means
df_means = behav.behavior_means(
    data_proc,
    by_conds=['modality', 'coherence', 'heading', 'delta'],
    long_format=True
)

# %% ================================================
# plot actual data
# ===================================================

# currently no fit curve, so will just draw lines between each
# can fit with a gaussian eventually (using behavior.utils.gauss_fit_hdg_group)

# behav.plot_behavior_hdg(
#     df_means[df_means['delta'] == 0],
#     col='coherence',
#     hue='modality',
#     palette=['k', 'r', 'b'],
#     hue_order=['ves', 'vis', 'comb'],
#     )

# %% ================================================
# Set up for model fitting
# ===================================================

data_delta0 = data[data['delta'] == 0]
X = data_delta0[["heading", "modality", "coherence", "delta"]]
y = data_delta0[["choice", "PDW", "RT"]]
y['RT'] += 0.3  # add offset for motion platform latency kluge

grid_vec = np.arange(-3, 0, 0.01) 
time_vec = np.arange(0, 2, 0.025)



init_params = {
    'kmult': [0.7, 0.7, 1.68],          # ves, vis sensitivites. if length 2, vis will be scaled by coh
    'bound': [1.0, 1.0, 0.95],          # ves, vis, comb bounds. 
    'non_dec_time': [0.13, 0.2, 0.5],   # non-decision time (secs)
    'wager_thr': [1.0, 1.0, 1.2],       # log odds threshold for high bets
    'wager_alpha': [0.03, 0.14, 0.05],  # base rate of low bets
    'cue_weights': [0.5, 0.5]
}

# initialize new DDM object
ddm_fit = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,  
    return_wager=True,
    )
SelfMotionDDM.save_params(init_params, save_dir / "init_params.json")

# set some options for the BADS routine (or scipy.minimize)
fit_options = {
    "random_seed": 0,
    "max_fun_evals": 500,
    "display": "full"
}

# ===================================================
# run the fit, with some fixed params
# ===================================================
ddm_fit.fit(
    X, y,
    fit_method='bads',
    fit_options=fit_options,
    )
ddm_fit.save(save_dir / "fitted_model.json")

trace_df = pd.DataFrame(ddm_fit.fit_trace_)
trace_df.to_csv(save_dir / "fit_trace.csv", index=False)

# print comparison table of params
df_params = pd.DataFrame([init_params, ddm_fit.params_])
df_params['name'] = ['Sim', 'Fit']
df_params.set_index('name', inplace=True)
print(df_params)

# save fitted params to disk
SelfMotionDDM.save_params(ddm_fit.params_, save_dir / "fitted_params.json")

# %% ================================================
# Plot fitted curves on top of original simulated data
# # ===================================================

# create a nicely spaced set of headings to run predictions over using fitted parameters
X_pred = dots3DMP_create_trial_list(
    hdgs=np.linspace(-12, 12, 100),
    mods=X["modality"].unique(),
    cohs=X["coherence"].unique(),
    nreps=1,            # only need 1
    shuff=False,        # shuffle is unnecessary here
)

# use the wager maps from fitting to the short list of headings,
# i.e. DON'T recompute them with the new spaced heading set!!
preds_, preds_samples = ddm_fit.predict(
    X_pred,
    n_samples=1,
    use_cached_wager_maps=True,
    rt_sampling_method="mean",
)

# a couple of kluges here to get the nice model predictions dataframe
# RT predictions come from the samples df, preds_ originally has the RT
# likelihoods needed for fitting
# for choice and PDW we can use preds_ columns as is, because they are
# probabilities of the binary outcome

# RT in task was defined as time from "motion onset",
# where motion onset is 1% of max acceleration
# in practice this is about 0.3 seconds after the stimulus onset
preds_['RT'] = preds_samples['RT'] - 0.3 
preds_full = process_predictions(X_pred, preds_)

actual_data_means = df_means[
    (df_means['delta'] == 0) & 
    (df_means['variable'].isin(['choice', 'PDW', 'RT']))
    ]
fitted_data = preds_full.loc[preds_full['delta'] == 0, :]

g = behav.plot_behavior_hdg(
    actual_data_means,
    data_fit=fitted_data,
    col='coherence',
    hue='modality',
    palette=['k', 'r', 'b'],
    hue_order=['ves', 'vis', 'comb'],
    )
g.figure.savefig(save_dir / "behavior_fit.png", dpi=150, bbox_inches="tight")
