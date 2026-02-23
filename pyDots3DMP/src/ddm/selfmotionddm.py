# %% ----------------------------------------------------------------

from copy import deepcopy
import itertools
import json
import logging
from pathlib import Path
import time
from typing import Any, Literal, Optional, Union

import numpy as np
import pandas as pd
from pybads import BADS
from scipy.optimize import minimize
from scipy.signal import convolve
from scipy.stats import norm, skewnorm, truncnorm

from behavior.utils import log_lik_bin, log_lik_cont, margconds_from_intersection

from .Accumulator import Accumulator

logger = logging.getLogger(__name__)

# %% ----------------------------------------------------------------
    
class SelfMotionDDM:
    PARAM_NAMES = ('kmult', 'bound', 'non_dec_time', 'wager_thr', 'wager_alpha')

    def __init__(
        self,
        grid_vec: np.ndarray,
        tvec: np.ndarray,
        kmult: list = [0.3, 0.3],
        bound: list = [1., 1., 1.],
        non_dec_time: list = [0.3],
        return_wager: bool = True,
        wager_thr: list = [1.],
        wager_alpha: list = [0.05],
        wager_maps: Optional[list] = None,
        wager_axis: Optional[int] = None,
        stim_scaling: Union[tuple[np.ndarray, np.ndarray], bool] = True,
        use_vectorized: bool = True,
        ):
        """3DMP accumulator model with confidence/wager readout
        :param grid_vec: vector of DV grid points
        :param tvec: vector of time points
        :param kmult: list of k multipliers [k_ves, k_vis]
        :param bound: list of bounds per modality [ves, vis, comb]      
        :param non_dec_time: list of non-decision times per modality
        :param return_wager: whether to compute/return wager predictions
        :param wager_thr: list of wager thresholds per modality
        :param wager_alpha: list of wager alpha parameters per modality
        :param wager_maps: list of existing wager maps to use for predictions
        :param wager_axis: axis for wager calculation (None = log odds)
        :param stim_scaling: whether to use stimulus-driven urgency signals
        :param use_vectorized: whether to use vectorized cdf/pdf implementations
        """
        self.grid_vec = grid_vec
        self.tvec = tvec
        self.kmult = kmult
        self.bound = bound
        self.non_dec_time = non_dec_time
        self.return_wager = return_wager
        self.wager_thr = wager_thr
        self.wager_alpha = wager_alpha
        self.wager_maps = wager_maps
        self.wager_axis = wager_axis
        self.stim_scaling = stim_scaling
        self.use_vectorized = use_vectorized

        # initialize internal containers used by fit/predict
        self.init_params = {k: getattr(self, k) for k in self.PARAM_NAMES}
        self.params_ = self.init_params.copy() # will hold all params after fitting
        self.fixed_params = {}

        self.accumulators_ = {}  # cache of accumulators per condition if desired

    def print_params(self, ndigits: int = 4) -> None:
        """Pretty-print params_ with rounded floats."""
        rounded = {k: round_val(v, ndigits) for k, v in self.params_.items()}
        for k, v in rounded.items():
            print(f"  {k}: {v}")

    def to_dict(self) -> dict:
        """Return a JSON-serializable dict of attributes needed to reinstantiate."""
        def to_serializable(v):
            if isinstance(v, np.ndarray):
                return v.tolist()
            if isinstance(v, (np.floating, np.integer)):
                return float(v) if isinstance(v, np.floating) else int(v)
            if isinstance(v, tuple) and len(v) == 2 and all(isinstance(x, np.ndarray) for x in v):
                return {"_tuple_arrays": [v[0].tolist(), v[1].tolist()]}
            if isinstance(v, list) and v and isinstance(v[0], np.ndarray):
                return [x.tolist() if isinstance(x, np.ndarray) else x for x in v]
            return v

        return {
            "grid_vec": to_serializable(self.grid_vec),
            "tvec": to_serializable(self.tvec),
            "params_": {k: to_serializable(v) for k, v in self.params_.items()},
            "return_wager": self.return_wager,
            "wager_maps": to_serializable(self.wager_maps) if self.wager_maps is not None else None,
            "wager_axis": self.wager_axis,
            "stim_scaling": to_serializable(self.stim_scaling) if isinstance(self.stim_scaling, tuple) else self.stim_scaling,
        }

    def save(self, path: Path | str) -> None:
        """Save the model state to a JSON file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: Path | str) -> "SelfMotionDDM":
        """Recreate a SelfMotionDDM from a JSON file saved by save()."""
        with open(path) as f:
            d = json.load(f)

        grid_vec = np.array(d["grid_vec"])
        tvec = np.array(d["tvec"])
        params_ = d["params_"]
        stim_scaling = d["stim_scaling"]
        if isinstance(stim_scaling, dict) and "_tuple_arrays" in stim_scaling:
            stim_scaling = tuple(np.array(a) for a in stim_scaling["_tuple_arrays"])
        wager_maps = d["wager_maps"]
        if wager_maps is not None:
            wager_maps = [np.array(w) for w in wager_maps]

        obj = cls(
            grid_vec=grid_vec,
            tvec=tvec,
            kmult=params_["kmult"],
            bound=params_["bound"],
            non_dec_time=params_["non_dec_time"],
            wager_thr=params_["wager_thr"],
            wager_alpha=params_["wager_alpha"],
            return_wager=d["return_wager"],
            wager_maps=wager_maps,
            wager_axis=d["wager_axis"],
            stim_scaling=stim_scaling,
        )
        obj.params_ = {k: v for k, v in params_.items()}
        return obj

    @classmethod
    def params_table(
        cls,
        *instances: "SelfMotionDDM",
        ndigits: int = 4,
    ) -> pd.DataFrame:
        """Build a DataFrame of params_ across one or more instances (rows=instances, columns=param names)."""
        if not instances:
            return pd.DataFrame(columns=cls.PARAM_NAMES)
        rows = [
            {p: round_val(inst.params_[p], ndigits) for p in cls.PARAM_NAMES}
            for inst in instances
        ]
        return pd.DataFrame(rows, columns=cls.PARAM_NAMES)

    def fit(
        self,
        X: pd.DataFrame,
        y: pd.DataFrame,
        fixed_params: Optional[list[str]]=None,
        fit_method: str = 'bads',
        fit_options: Optional[dict] = None,
        ) -> 'SelfMotionDDM':
        """fit model to data in X and y, with optional fixed parameters"""

        logger.info('Starting model fitting')

        fit_start_time = time.perf_counter()
        
        self.n_features_in_ = len(X.columns)

        # get list of fittable parameters, if any
        self.params_ = self.init_params.copy() 
        params_list = self._get_fit_params(fixed_params)

        if params_list:
            # concatenate into single array for passing to optimization function
            # but store original list lengths for reconstructing dict later
            params_array = np.array(list(itertools.chain(*params_list)))
            self.param_end_inds = list(itertools.accumulate(map(len, params_list)))

            # pass data as fixed inputs to objective function
            optim_fcn_part = lambda params: self._objective_fcn(params, X, y)

            if self.save_dir is not None:
                Path.mkdir(self.save_dir, exist_ok=True)
            if fit_method.lower() == 'bads':

                # TODO expose these to the user
                lb = params_array * 0.25
                ub = params_array * 3.0
                plb = params_array * 0.5
                pub = params_array * 2.0
                bads_bounds = (lb, ub, plb, pub)
                
                bads = BADS(
                    optim_fcn_part, 
                    params_array, 
                    *bads_bounds, 
                    options=fit_options
                    )
                result = bads.optimize()

            else:
                result = minimize(
                    self._objective_fcn,
                    params_array,
                    args=(X, y),
                    method=fit_method,
                    options=fit_options
                    )

            logger.info("================")
            logger.info(result)
            logger.info("================")
                   
            # at the end, store fitted params back into dict, with fixed ones
            self._build_params_dict(result.x, self.param_end_inds)

        fit_duration = time.perf_counter() - fit_start_time
        logger.info(f"Fitting took {fit_duration:3f}s / {fit_duration/60:3f} mins")
        
        return self

    def _objective_fcn(
        self,
        params_array: np.ndarray,
        X: pd.DataFrame,
        y: pd.DataFrame,
        fixed_params: Optional[dict[str, Any]]=None,
        ) -> float:
        """
        objective function for optimization - negative log likelihood of data given model params
        return float for minimization
        """

        params_array = np.asarray(params_array)
            
        # combine params array passed to objective function with fixed params
        # to reconstruct full params dict expected by custom predict method
        self._build_params_dict(params_array, self.param_end_inds, fixed_params)
        # logger.info(params_array)
        # print('Current params: %s', {k: [round(vv, 2) for vv in v] for k, v in self.params_.items()})

        t0_pred = time.perf_counter()
        y_pred, _ = self.predict(X, y)
        t1_pred = time.perf_counter() - t0_pred
        logger.debug(f'single objective function prediction run took {t1_pred:.2f} seconds')
        
        # calculate log likelihoods for each output
        log_lik_choice = log_lik_bin(y['choice'].to_numpy(), y_pred['choice'].to_numpy()) / len(y)
        log_lik_pdw    = log_lik_bin(y['PDW'].to_numpy(), y_pred['PDW'].to_numpy()) / len(y)
        log_lik_rt     = log_lik_cont(y_pred['RT'].to_numpy()) / len(y)
        
        self.log_lik_ = {
            'choice': log_lik_choice,
            'pdw': log_lik_pdw,
            'rt': log_lik_rt
        }
        if self.return_wager:
            logger.debug('Log likelihoods - choice: %.2f, PDW: %.2f, RT: %.2f', 
                        log_lik_choice, log_lik_pdw, log_lik_rt)
            self.neg_llh_ = -sum([log_lik_choice, log_lik_pdw, log_lik_rt])
        else:
            logger.debug('Log likelihoods - choice: %.2f, RT: %.2f', 
                        log_lik_choice, log_lik_rt)
            self.neg_llh_ = -sum([log_lik_choice, log_lik_rt])
        logger.debug('Total loss:\t%.2f', self.neg_llh_)

        return self.neg_llh_

    def predict(
        self,
        X,
        y=None,
        n_samples: int = 1,
        cache_accumulators: bool = False,
        use_cached_wager_maps: bool = False,
        rt_sampling_method: Literal["sample", "mean", "mode"] = "sample",
        seed=None
        ):
        """
        generate model predictions for data in X
        :param X: DataFrame with columns modality, coherence, delta, heading
        :param y: DataFrame with columns choice, PDW, RT (for RT likelihood calculation)
        :param n_samples: number of samples to draw for probabilistic predictions (0 = none)
        :param cache_accumulators: (default = False) whether to cache accumulator objects
        :param use_cached_wager_maps: (default = True) whether to use existing cached_wager_maps
        :param rt_sampling_method: (default="sample) how to draw predicted RTs from distribution - options are "sample", "mean", or "mode"
        :param seed: random seed for sampling
        :return: 
            predictions - DataFrame with columns choice, PDW, RT (predicted likelihoods)
            pred_sample - DataFrame with sampled predictions, returns None if n_samples == 0 
                (n_samples > 0 -> 
                    how many draws from binomial for choice/PDW, or how many samples from RT_dist for RT
                    if rt_sampling_method is "mean" or "mode", will instead take expected value or argmax of RT_dist
                )
        """
        
        rng = np.random.RandomState(seed)
        
        mods = np.unique(X['modality']).astype(float)
        cohs = np.unique(X['coherence']).astype(float)
        deltas = np.unique(X['delta']).astype(float)
        hdgs, hdg_inds = np.unique(X['heading'], return_inverse=True)
        hdgs = hdgs.astype(float)

        K_SCALE_FACTOR = 1e3
        if not self.stim_scaling:
            b_ves, b_vis = np.ones_like(self.tvec), np.ones_like(self.tvec)
        elif isinstance(self.stim_scaling, tuple):
            b_ves, b_vis = self.stim_scaling
        else:
            b_ves, b_vis = get_stim_urgs(self.tvec)
            K_SCALE_FACTOR = 1e4

        b_ves /= len(b_ves)
        b_vis /= len(b_vis)
        b_vals = [b_ves, b_vis, np.vstack((b_ves, b_vis)).T]

        # handle parameters per modality
        kves, kvis = self._handle_kmult(self.params_['kmult'], cohs.T, k_scale=K_SCALE_FACTOR) 
        bound = self._handle_param_mod(self.params_['bound'], mods)  
        non_dec_time = self._handle_param_mod(self.params_['non_dec_time'], mods)  
        thetas = self._handle_param_mod(self.params_['wager_thr'], mods)  
        alphas = self._handle_param_mod(self.params_['wager_alpha'], mods)  

        # initialize predictions dataframes
        predictions = pd.DataFrame(
            np.full((X.shape[0], 3), fill_value=np.nan), columns=['choice', 'PDW', 'RT']
            )
        pred_sample = deepcopy(predictions) if n_samples else None

        # compute wager maps if not existing, or cache not requested
        if self.return_wager and (not use_cached_wager_maps or not self.wager_maps):

            self.wager_maps = []

            # ves, vis, comb overall sensitivities
            k_vals_fixed = [kves, kvis.mean().item(), [kves, kvis.mean().item()]]
            
            for m, mod in enumerate(mods):

                # set accumulators with absolute drifts, for log odds mappings
                accumulator = Accumulator(grid_vec=self.grid_vec, tvec=self.tvec, bound=bound[m])

                abs_drifts, t_eff = calc_selfmotion_drifts(
                    b_vals[m], k_vals_fixed[m], self.tvec, hdgs[hdgs>=0], delta=0,
                    )
                accumulator.tvec = t_eff

                # run the method of images - diffusion to bound to extract pdfs, cdfs, and LPO
                accumulator.apply_drifts(abs_drifts, hdgs[hdgs>=0])
                accumulator.compute_distrs(return_pdf=True, use_vectorized=self.use_vectorized) # get the pdfs for wager calculation

                if cache_accumulators:
                    self.accumulators_[('wager', mod)] = accumulator

                if self.wager_axis is None:
                    log_odds_map = accumulator.log_posterior_odds()
                    self.wager_maps.append(log_odds_map)
                else:
                    raise NotImplementedError('alternatives to log odds not yet implemented')

        # boolean mask on wager map for high bets
        wager_is_high = [p >= theta for p, theta in zip(self.wager_maps, thetas)]

        # now loop over coherences and deltas with one accumulator each for actual predictions
        for c, coh in enumerate(cohs):
            k_vals = [kves, kvis[c], [kves, kvis[c]]]

            for m, mod in enumerate(mods):
                for d, delta in enumerate(deltas):

                    trial_index = (X['modality'] == mod) & (X['coherence'] == coh) & (X['delta'] == delta)
                    trial_index = trial_index.to_numpy()

                    # non-valid conditions
                    if (delta != 0 and mod < 3) or (c > 0 and mod == 1) or trial_index.sum()==0:
                        continue       

                    # set up accumulator for this condition
                    accumulator = Accumulator(grid_vec=self.grid_vec, tvec=self.tvec, bound=bound[m])
                    drifts, t_eff = calc_selfmotion_drifts(
                        b_vals[m], k_vals[m], self.tvec, hdgs, delta=delta,
                        )
                    accumulator.tvec = t_eff
                    
                    # this time use signed headings
                    accumulator.apply_drifts(drifts, hdgs) 

                    # run the method of images - diffusion to bound to extract pdfs, cdfs, and LPO
                    accumulator.compute_distrs(return_pdf=self.return_wager, use_vectorized=self.use_vectorized)

                    if cache_accumulators:
                        self.accumulators_[(mod, coh, delta)] = accumulator
                    # get predictions for all trials in this condition
                    
                    # ====== CHOICE ======
                    p_right = np.clip(accumulator.p_corr_.T, 1e-10, 1-1e-10)
                    predictions.loc[trial_index, 'choice'] = p_right[hdg_inds[trial_index]]

                    # ====== WAGER ======
                
                    for h, hdg in enumerate(hdgs):
                        
                        trial_index = (X['modality'] == mod) & (X['coherence'] == coh) & \
                                    (X['heading'] == hdg) & (X['delta'] == delta)
                        trial_index = trial_index.to_numpy()
                        if trial_index.sum() == 0:
                            continue

                        p_choice = np.array([p_right[h], 1 - p_right[h]])
                        
                        # # ====== CHOICE ======                        
                        if n_samples:
                            pred_sample.loc[trial_index, 'choice'] = rng.binomial(n_samples, p_right[h], trial_index.sum()) / n_samples

                        # ====== WAGER ======
                        if self.return_wager:
                            # select pdf for losing race, given correct or incorrect
                            pxt_up = np.squeeze(accumulator.up_lose_pdf_[h, :, :])
                            pxt_lo = np.squeeze(accumulator.lo_lose_pdf_[h, :, :])
                            total_p = np.sum(pxt_up + pxt_lo)
                            pxt_up /= total_p
                            pxt_lo /= total_p

                            p_choice_and_wager = np.array(
                                [
                                    [
                                        np.sum(pxt_up[wager_is_high[m]]),   # pRight+High
                                        np.sum(pxt_up[~wager_is_high[m]])   # pRight+Low
                                    ],   
                                    [
                                        np.sum(pxt_lo[wager_is_high[m]]),   # pLeft+High
                                        np.sum(pxt_lo[~wager_is_high[m]])   # pLeft+Low
                                    ]
                                ]
                            )

                            # calculate p_wager using Bayes rule, then factor in base rate of low bets ("alpha")
                            p_choice_given_wager, p_wager = margconds_from_intersection(
                                p_choice_and_wager, p_choice
                                )
                            p_wager += np.array([-alphas[m], alphas[m]]) * p_wager[0]
                            p_wager = np.clip(p_wager, 1e-100, 1-1e-100)

                            predictions.loc[trial_index, 'PDW'] = p_wager[0] # proportion of high bets

                            if n_samples:
                                pred_sample.loc[trial_index, 'PDW'] = rng.binomial(n_samples, p_wager[0], trial_index.sum()) / n_samples


                        # ====== RT ======

                        # first convolve model RT distribution with non-decision time
                        ndt_dist = norm.pdf(self.tvec, loc=non_dec_time[m], scale=0.2) #scale=self.params_['sigma_ndt'])
                        rt_dist = np.squeeze(accumulator.rt_dist_[h, :])
                        rt_dist = convolve(rt_dist, ndt_dist / ndt_dist.sum())
                        rt_dist = np.clip(rt_dist, 1e-10, a_max=None)

                        # trim to original length of tvec and renormalize to get posterior
                        rt_dist = rt_dist[:len(accumulator.tvec)]
                        rt_dist /= rt_dist.sum()

                        if y is not None:
                            # need original data here to get RT likelihoods
                            actual_rts = y.loc[trial_index, 'RT'].values
                            # dist_inds = [np.argmin(np.abs(self.tvec - rt)) for rt in actual_rts]
                            dist_inds = np.searchsorted(self.tvec, actual_rts)
                            dist_inds[dist_inds >= len(self.tvec)] = len(self.tvec) - 1  # cap at max index
                            predictions.loc[trial_index, 'RT'] = rt_dist[dist_inds]

                        if n_samples:
                            if rt_sampling_method == "sample":
                                sampled_RTs = np.random.choice(
                                    self.tvec, (trial_index.sum(), n_samples), replace=True, p=rt_dist
                                    )
                                pred_sample.loc[trial_index, 'RT'] = sampled_RTs.mean(axis=1)
                            elif rt_sampling_method == "mean":
                                pred_sample.loc[trial_index, 'RT'] = np.dot(self.tvec, rt_dist)
                            elif rt_sampling_method == "mode":
                                pred_sample.loc[trial_index, 'RT'] = self.tvec[np.argmax(rt_dist)]

        return predictions, pred_sample
    

    def _build_params_dict(
        self,
        params_array: np.ndarray,
        end_inds: list,
        fixed_params: Optional[dict[str, Any]]=None,
        ) -> None:
        """reconstruct full params dictionary from params array and fixed params dictionary"""

        fixed_params = fixed_params if fixed_params is not None else self.fixed_params

        start_ind = 0
        for p, param in enumerate(self.fit_param_names):
            if p > 0:
                start_ind = end_inds[p-1]
            self.params_[param] = params_array[start_ind:end_inds[p]].tolist()
        self.params_ = self.params_ | fixed_params

        # for pn, p in self.params_.items():
        #     print(f'{pn}:', sep='\t')
        #     for elem in p:
        #         print(f'{elem:.2f}', sep=',')


    def _get_fit_params(self, fixed_params) -> list:
        """return list of fittable parameter values (from initial parameters dict)"""

        if fixed_params:
            self.fixed_params = {k: self.init_params[k] for k in fixed_params}
            
        self.fit_param_names = [k for k in self.init_params.keys() if k not in self.fixed_params.keys()]
        logger.debug(f'Fixed parameters: {self.fixed_params}')
        logger.debug(f'Fitting parameters: {self.fit_param_names}')

        return [self.init_params[k] for k in self.fit_param_names]


    def simulate(
        self,
        X,
        n_samples: int = 1,
        sample_dvs: bool = False,
        seed=None):
        """
        simulate data from model given X
        :param X: DataFrame with columns modality, coherence, delta, heading
        :param n_samples: number of samples to draw per trial (default 1)
        :param return_dvs: whether to sample DV trajectories
        :param seed: random seed for sampling
        :return: simulated DataFrame with columns choice, PDW, RT

        """
        
        _, preds = self.predict(
            X,
            y=None,
            n_samples=n_samples, 
            cache_accumulators=sample_dvs, 
            seed=seed
            )

        if sample_dvs:
            # build n_samples DV trajectories for unique conditions (from cached accumulator objects)

            # NOTE: IMCOMPLETE
            
            Xu = X.drop_duplicates().values
            preds = pd.DataFrame(np.repeat(Xu, n_samples, axis=0), columns=X.columns)

            mods = np.unique(X['modality']).astype(float)
            non_dec_time = self._handle_param_mod(self.params_['non_dec_time'], mods)  
            alphas = self._handle_param_mod(self.params_['wager_alpha'], mods)  
            
            for i, row in X.iterrows():
                mod = int(row['modality'])
                coh = row['coherence']
                delta = row['delta']

                accum = self.accumulators_[(mod, coh, delta)]

                ndt_mean = non_dec_time[mod-1]
                ndt_std = 0.050
                ndt_min = ndt_mean / 2
                ndt_max = ndt_mean + ndt_min
                
                for ihdg, hdg in enumerate(accum.drift_labels):

                    trial_inds = preds.index[
                        (preds["modality"]==mod) & (preds["coherence"]==coh) & \
                            (preds["delta"]==delta) & (preds["heading"]==hdg)
                    ].to_numpy()
                    
                    ndt = truncnorm.rvs(ndt_min, ndt_max, loc=ndt_mean, scale=ndt_std, 
                                    size=len(trial_inds))

                    for itr in range(n_samples):
                        dv = accum.dv(ihdg)

                        is_hit_bnd = (dv >= accum.bound).any(axis=0)
                        t_bnd_cross = np.argmax((dv >= accum.bound) == 1, axis=0)

                        if is_hit_bnd.all():
                            # both accumulators hit the bound - choice and RT determined by whichever hit first
                            choice = np.argmin(t_bnd_cross)
                            rt_ind = t_bnd_cross[choice]
                            final_v = dv[rt_ind, choice ^ 1]  # losing accumulator at RT

                        elif ~is_hit_bnd.any():
                            # neither accumulator hits the bound
                            rt_ind = -1
                            choice = np.argmax(np.argmax(dv, axis=0))  # winner is whoever has max dv value
                            final_v = accum.bound[choice] - (dv[rt_ind, choice] - dv[rt_ind, choice ^ 1])
                            # wager odds map accounts for the distance between winner and loser, so this shifts up both accumulators
                            # as if the 'winner' did hit the bound, so we can do a consistent look-up on the wager odds map

                        else:
                            # only one hits the bound
                            choice = np.argmax(is_hit_bnd)
                            rt_ind = t_bnd_cross[choice]
                            final_v = dv[rt_ind, choice ^ 1]

                        if self.return_wager:

                            # look-up log odds threshold
                            grid_ind = np.argmin(np.abs(accum.grid_vec - final_v))

                            wager_accum = self.accumulators_[('wager', mod)]
                            
                            # log_odds = wager_odds_maps[m][rt_ind, grid_ind]
                            wager = int(wager_accum.wager_map[rt_ind, grid_ind])
                            wager *= (np.random.random() > alphas[mod-1])  # incorporate base-rate of low bets
                            preds.loc[trial_inds[itr].item(), 'PDW'] = wager

                        # flip choice result so that left choices = 0, right choices = 1 in the output
                        preds.loc[trial_inds[itr].item(), 'choice'] = choice ^ 1

                        # RT = decision time + non-decision time
                        preds.loc[trial_inds[itr].item(), 'RT'] = self.tvec[rt_ind] + ndt[itr]

        return preds
        
    
    @staticmethod
    def _handle_kmult(kmult, cohs, k_scale=1):
        """break down kmult into (kves, kvis)"""
        kmult = [k*k_scale for k in kmult]
        kves = kmult[0]
        if len(kmult) == 1:
            kvis = kmult * cohs
            kves = np.mean(kvis)  # ves lies between vis cohs
        if len(kmult) == 2:
            kvis = kmult[1] * cohs # vis scaled by coherence
        else:
            kvis = kmult[1:]  # 3 independent kmults
        return kves, kvis

    @staticmethod
    def _handle_param_mod(param, mods):
        """broadcast param to match number of mods"""
        n_mods = len(mods)
        
        if np.isscalar(param):
            return [param] * n_mods

        param_list = list(param)
        if len(param_list) == n_mods:
            return param_list
        if len(param_list) == 1:
            return param_list * n_mods

        raise ValueError("Length of param list does not match number of modalities")

# %% -----------------------

def get_stim_urgs(
    tvec: np.ndarray,
    pos: Optional[np.ndarray] = None,
    skew_params: Optional[tuple] = None
    ):
    """
    Return acceleration and velocity profiles for stimulus weighting

    Args:
        tvec (np.ndarray): time vector, used to create position profile with skew_params. 
                            Ignored if pos provided directly.
        pos (Optional[np.ndarray], optional): position profile. Defaults to None.
        skew_params (Optional[tuple], optional): args for skewnorm.cdf. Defaults to (2, 0.8, 0.4).

    Returns:
        np.ndarray, np.ndarray: acceleration, velocity vectors
    """
    if pos is None:
        ampl = 0.16

        # pos = norm.cdf(tvec, 0.9, 0.3) * ampl
        if skew_params is None:
            skew_params = (2, 0.8, 0.4)
        pos = skewnorm.cdf(tvec, *skew_params) * ampl

    vel = np.gradient(pos)
    acc = np.gradient(vel)

    vel /= vel.max()
    acc /= acc.max()
    
    return acc, vel

    
def calc_selfmotion_drifts(
    b_t: np.ndarray,
    b_k: Union[float, tuple[float, float]],
    tvec: float,
    hdgs: np.ndarray,
    delta: float = 0.0, 
    cue_weights: Optional[tuple[float, float]] = None
    ) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate instaneous drift rates given time course and modality sensitivities
    See Drugowitsch et al 2014 eLife supplementary equations

    b_t:  time-course sensitivity
    b_k:  stimulus modality sensitivity (tuple --> bimodal condition)
    tvec: original (linear) time vector
    hdgs: heading values
    delta: heading delta (combined modality)
    cue_weights: tuple for custom override of cue weights
        (default is None, in which case optimal cue weights are computed from b_k)

    Returns:
       drifts - array of instantaneous drift rates (T x drifts)
       t_eff - effective time (momentary power of stimulus)
       
    """

    sin_hdgs = np.sin(np.deg2rad(hdgs))
    dt = np.gradient(tvec)

    cumul_bt = np.cumsum((b_t**2)/(b_t**2).sum(axis=0), axis=0)
    
    if isinstance(b_k, (int, float)):
        # only one sensitivity and time-course - ves or vis (logic is the same)
        t_eff = cumul_bt * tvec[-1]
        drifts = np.reshape(b_t, (-1, 1))**2 * b_k * sin_hdgs # see Drugowitsch et al. 2014 supp eq 2 & 7
        # drifts = b_k * sin_uhdgs   # w/o stim scaling, reduces to this

    elif len(b_k) == 2:
        # two sensitivities/time-courses - combined condition

        b_k = np.array(b_k, dtype=float)
        if cue_weights is None:
            k2 = b_k**2
            cue_weights = np.sqrt(k2 / k2.sum())
            
        w_ves, w_vis = cue_weights
        
        # +ve delta means ves to the left, vis to the right
        drift_ves = b_t[:, [0]]**2 * b_k[0] * np.sin(np.deg2rad(hdgs - delta / 2))
        drift_vis = b_t[:, [1]]**2 * b_k[1] * np.sin(np.deg2rad(hdgs + delta / 2))
        
        drifts = w_ves * drift_ves + w_vis * drift_vis

        # Drugowitsch et al. 2014 supp eq 14
        t_eff = (w_ves**2 * cumul_bt[:,0] + w_vis**2 * cumul_bt[:,1]) * tvec[-1]

    # cumsum over time, then divide by t_eff to get instantaneous drifts
    drifts = np.cumsum(drifts, axis=0) / t_eff[:, None]

    return drifts, t_eff


# param printing util
def round_val(v, ndigits):
    if isinstance(v, np.ndarray):
        v = v.tolist()
    if isinstance(v, (list, tuple)):
        return [round(float(x), ndigits) if isinstance(x, (int, float, np.floating)) else x for x in v]
    if isinstance(v, (int, float, np.floating)):
        return round(float(v), ndigits)
    return v