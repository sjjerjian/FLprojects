"""Fit some real data"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

# custom imports
from ddm import SelfMotionDDM
from behavior.utils import dots3DMP_create_trial_list, data_cleanup
import behavior.descriptive as behav

from datetime import datetime

# %% ===== set up save location and logger =====
save_dir = Path(f"results/lucio/{datetime.now().strftime('%y%m%d_%H%M%S')}")
Path.mkdir(save_dir, parents=True, exist_ok=True)

# Set up logger to console and file - this will show us the logging info for our DDM, and the bads fitting routine
logger = logging.getLogger()
logger.setLevel(logging.INFO)
fmt = logging.Formatter("%(asctime)s | %(name)s | %(levelname)s | %(message)s", datefmt="%y%m%d_%H%M%S")
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    
    fh = logging.FileHandler(save_dir / "fit_lucio.log", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)


# %% ===============
# load data 
# ==================

datafilepath = "/Users/stevenjerjian/FLprojects/lucio_20220512-20230606.csv"
data = data_cleanup(datafilepath)  # this is a hack function to quickly load and clean the data, could be improved/generalized with options

# %% set up for fitting

def process_predictions(X, y=None):

    mod_map = {1: "ves", 2: "vis", 3: "comb"}
    data = pd.concat((X, y), axis=1) if y is not None else X.copy()

    # replicate ves for high coherence, for plotting convenience
    data = behav.replicate_ves(data) 
    # map modalities from ordinal to str labels
    data['modality'] = data['modality'].map(mod_map)

    return data

data_proc = process_predictions(data)
# plot_behavior_hdg expects pre-computed means
df_means = behav.behavior_means(
    data_proc,
    by_conds=['modality', 'coherence', 'heading', 'delta'],
    long_format=True)

# %%
# plot simulated data. currently no fit curve, so will just draw lines between each
# can fit with a gaussian eventually (using behavior.utils.gauss_fit_hdg_group)
behav.plot_behavior_hdg(
    df_means[df_means['delta'] == 0],
    col='coherence',
    hue='modality',
    palette=['k', 'r', 'b'],
    hue_order=['ves', 'vis', 'comb'],
    )

# %% RUN MODEL FITTING

data_delta0 = data[data['delta'] == 0]
X = data_delta0[["heading", "modality", "coherence", "delta"]]
y = data_delta0[["choice", "PDW", "RT"]]
y['RT'] += 0.3  # add offset for motion platform latency kluge

grid_vec = np.arange(-3, 0, 0.01) 
time_vec = np.arange(0, 2, 0.025)

init_params = {
    'kmult': [1.5, 2.2],          # ves, vis sensitivites. if length 2, vis will be scaled by coh
    'bound': [1.0],               # ves, vis, comb bounds. 
    'non_dec_time': [0.1],      # non-decision time (secs)
    'wager_thr': [1.6],         # log odds threshold for high bets
    'wager_alpha': [0.06],      # base rate of low bets
}
init_params = {
    'kmult': [0.7530327, 1.55471721, 2.80188416],          # ves, vis sensitivites. if length 2, vis will be scaled by coh
    'bound': [0.96068757, 0.51813231, 0.86937995],               # ves, vis, comb bounds. 
    'non_dec_time': [0.07577165, 0.37180706, 0.05],      # non-decision time (secs)
    'wager_thr': [0.93520996, 1.03769139, 1.06222401],         # log odds threshold for high bets
    'wager_alpha': [0.02533596, 0.10985811, 0.09932743],      # base rate of low bets
}

with open(save_dir / "init_params.json", "w") as f:
    json.dump(init_params, f, indent=4)
    
# initialize new DDM object
ddm_fit = SelfMotionDDM(
    grid_vec=grid_vec,
    tvec=time_vec,
    **init_params, 
    stim_scaling=True,  
    return_wager=True,
    )

# ddm_fit = SelfMotionDDM.load("fit_lucio_20260423_223747/fitted_model_lucio_freemods.json")
# print(SelfMotionDDM.params_table(ddm_fit))

# set some options for the BADS routine (or scipy.minimize)
fit_options = {
        "random_seed": 42,
        "max_fun_evals": 500,
        "display": "full"
    }

# %% run the fit, with some fixed params
ddm_fit.fit(
    X, y,
    fit_method='bads',
    fit_options=fit_options
    )
ddm_fit.save(save_dir / "fitted_model.json")

trace_df = pd.DataFrame(ddm_fit.fit_trace_)
trace_df.to_csv(save_dir / "fit_trace.csv", index=False)

# print comparison table of params
df_params = pd.DataFrame([init_params, ddm_fit.params_])
df_params['name'] = ['Sim', 'Fit']
df_params.set_index('name', inplace=True)
print(df_params)

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

# RT in task was defined as time from "motion onset" where motion onset is 1% of max acceleration
# in practice this is about 0.3 seconds after the stimulus onset
preds_['RT'] = preds_samples['RT'] - 0.3 
preds_full = process_predictions(X_pred, preds_)
# %%
g = behav.plot_behavior_hdg(
    df_means[(df_means['delta'] == 0) & (df_means['variable'].isin(['choice', 'PDW', 'RT']))],
    data_fit=preds_full.loc[preds_full['delta'] == 0, :],    # model fit for delta = 0
    col='coherence',
    hue='modality',
    palette=['k', 'r', 'b'],
    hue_order=['ves', 'vis', 'comb'],
    )
g.figure.savefig(save_dir / "behavior_fit.png", dpi=150, bbox_inches="tight")
# %%
