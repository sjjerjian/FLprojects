# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import seaborn as sns

# from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit

from functools import wraps
from typing import Optional
from .utils import prop_se, cont_se, gaus


# TODO functions list
# 1. extracting behavior from neural recording dataset
# 2. RT quantiles, plot vs confidence/accuracy
# 3. correct vs error metrics
# 4. decorator to loop any function over a grouping variable (i.e. day/block)
# 5. regression analyses


def behavior_means(df, by_conds='heading', drop_na=True, long_format=True):

    """
    Compute the mean and standard error of the data for each condition.
    Args:
        df (pd.DataFrame): the dataframe to compute the means of
        by_conds (list): the conditions to group by
        drop_na (bool): whether to drop rows with missing values
        long_format (bool): whether to return the data in long format
    Returns:
        pd.DataFrame: the dataframe with the means and standard errors
    """
    
    agg_funcs = {
        'choice': ['count', 'mean', prop_se],
        'PDW': ['count', 'mean', prop_se],
        'RT': ['count', 'mean', cont_se],
        'correct': ['count', 'mean', prop_se],
    }
    agg_funcs = {k:v for k,v in agg_funcs.items() if k in df.columns}

    output_vars = list(agg_funcs.keys())

    df_means = df.groupby(by=by_conds)[output_vars].agg(agg_funcs)
    if drop_na:
        df_means.dropna(axis=0, inplace=True)
    df_means = df_means.reset_index()
    df_means.columns = ['_'.join(col) if col[0] in output_vars else col[0] for col in df_means.columns]  # remove multi-level index

    if long_format:
        count_melt = df_means.melt(
            id_vars=by_conds,
            value_vars=df_means.columns[df_means.columns.str.contains('count')],
            var_name='variable', value_name='count')

        means_melt = df_means.melt(
            id_vars=by_conds,
            value_vars=df_means.columns[df_means.columns.str.contains('mean')],
            var_name='variable', value_name='mean')

        sems_melt = df_means.melt(
            id_vars=by_conds,
            value_vars=df_means.columns[df_means.columns.str.contains('_se')],
            var_name='variable', value_name='se')

        df_means = count_melt.copy()
        df_means['mean'] = means_melt['mean']
        df_means['se'] = sems_melt['se']

        df_means['variable'] = df_means['variable'].apply(lambda x: x.split('_')[0])

    return df_means


def replicate_ves(df):
    """Replicate vestibular condition rows, for every coherence level."""
    if 'coherence' in df.columns and len(np.unique(df['coherence'])) > 1:
        ucohs = np.unique(df['coherence'])

        dup_ves = df.loc[df['modality'] == 1, :]
        ucohs = np.unique(df['coherence'])
        ucohs = ucohs[ucohs != np.unique(dup_ves['coherence'])]

        result_df = pd.concat([dup_ves.assign(coherence=coh) for coh in ucohs],
                              ignore_index=True)
        result_df = pd.concat((df, result_df), ignore_index=True)
    else:
        result_df = df.copy()

    return result_df


def logit_fit_choice_hdg(df, num_hdgs: int = 200) -> pd.Series:

    hdgs = np.unique(df['heading']).reshape(-1, 1)
    xhdgs = np.linspace(np.min(hdgs), np.max(hdgs), num_hdgs).reshape(-1, 1)

    # logreg = smf.logit("choice ~ heading", data=df).fit()  # formula api
    logreg = sm.Logit(df['choice'], sm.add_constant(df['heading'])).fit()
    yhat = logreg.predict(sm.add_constant(xhdgs))
    params = logreg.params['heading'], logreg.params['const']

    # alternatively, using sklearn # TODO test
    # logreg = LogisticRegression().fit(df[['heading]], df['choice'])
    # yhat = logreg.predict_proba(xhdgs)[:, 1]

    return pd.Series({'yhat': yhat, 'params': params})


def gauss_fit_hdg(df, p0: np.ndarray, y_var: str = 'choice', numhdgs: int = 200) -> pd.Series:
    """Gaussian fitting over headings to single set of data."""
    hdgs = np.unique(df['heading']).reshape(-1, 1)
    xhdgs = np.linspace(np.min(hdgs), np.max(hdgs), numhdgs).reshape(-1, 1)

    if y_var == 'choice':
        probreg = sm.Probit(df['choice'], sm.add_constant(df['heading'])).fit(start_params=p0)
        yhat = probreg.predict(sm.add_constant(xhdgs))
        params = probreg.params['heading'], probreg.params['const']

    elif y_var == 'PDW':
        params, pcov = curve_fit(gaus, xdata=df['heading'], ydata=1-df['PDW'], p0=p0)
        yhat = 1 - gaus(xhdgs, *params).flatten()

    elif y_var == 'RT':
        params, pcov = curve_fit(gaus, xdata=df['heading'], ydata=df['RT'], p0=p0)
        yhat = gaus(xhdgs, *params).flatten()

    elif y_var == 'correct':
        raise NotImplementedError("Gaussian fit to p(correct) not implemented yet")

    # # To compute 1SD error on parameters,
    # perr = np.sqrt(np.diag(pcov))

    return pd.Series({'hdgs': xhdgs.flatten(), 'yhat': yhat, 'params': params})


def gauss_fit_hdg_group(
    df: pd.DataFrame, 
    p0: np.ndarray, 
    y_vars: tuple = ('choice', 'PDW', 'RT'), 
    by_conds: str = 'modality', 
    numhdgs: int = 200
    ) -> dict:
    """
    Fit a Gaussian to the data for each condition.
    Args:
        df (pd.DataFrame): the dataframe to fit the Gaussian to
        p0 (np.ndarray): the initial parameters for the Gaussian
        y_vars (tuple): the variables to fit the Gaussian to
        by_conds (str): the condition to group by
        numhdgs (int): the number of headings to fit the Gaussian to
    Returns:
        dict: a dictionary of the fit results
    """
    
    fit_results = {
        y_var: df.groupby(by=by_conds).apply(gauss_fit_hdg, p, y_var, numhdgs).dropna(axis=0).reset_index()
        for y_var, p in zip(y_vars, p0)
    }

    # explode hdgs and yhat vectors
    fit_results = {k: df.drop('params', axis=1).explode(['hdgs', 'yhat']) for k, df in fit_results.items()}

    return fit_results


def fit_results_to_dataframe(fit_results, by_conds = 'modality'):

    y_vars = list(fit_results.keys())

    by_conds.append('hdgs')
    fit_df = fit_results[y_vars[0]][by_conds].rename(columns={'hdgs': 'heading'})

    for label, df in fit_results.items():
        fit_df[label] = df['yhat']

    return fit_df


# %%

def cue_weighting(fit_results):

    wves_emp = None
    wves_pred = None

    for res in fit_results.keys():
        ...
        # TODO calculate wves pred and wves emp for each of pRight, PDW, RT
    return wves_emp, wves_pred

# %%


def plot_behavior_hdg(
    data_obs,
    data_fit: Optional[pd.DataFrame] = None,
    row: str = 'variable',
    col: str ='coherence',
    hue: str = 'modality',
    palette = sns.color_palette(),     
    hue_order: Optional[list] = None,
    **fig_kwargs
    ):  

    def _errbar_plot(x, y, yerr, **kwargs):
        plt.errorbar(x, y, yerr, **kwargs)

    # plot the empirical data points, using FacetGrid for convenient conditional plotting
    # see https://seaborn.pydata.org/generated/seaborn.FacetGrid.html
    g = sns.FacetGrid(
        data_obs,
        row=row,
        col=col,
        hue=hue,
        hue_order=hue_order,
        palette=palette,
        sharey=False,
        **fig_kwargs
        )   
    
    line_style = '-' if data_fit is None else ''
    g.map_dataframe(
        _errbar_plot, 'heading', 'mean', 'se', linestyle=line_style, marker='.'
        )

    # Single legend: three colors, same order as curves and points
    legend_handles = [
        Line2D([0], [0], color=c, lw=2, label=name)
        for c, name in zip(palette, hue_order)
    ]
    
    # overlay the fit data as a line
    for iax, (ax_key, ax) in enumerate(g.axes_dict.items()):
        if data_fit is not None and col is not None:
            ax_data = data_fit.loc[data_fit[col]==ax_key[1], :]
            y = ax_key[0]

            sns.lineplot(
                data=ax_data,
                x='heading',
                y=y,
                hue=hue,
                hue_order=hue_order,
                ax=ax,
                palette=palette,
                legend=False
                )

        ax.set_title("")
        ax.set_xlabel("")
        if 'choice' in ax_key:
            ax.set_title(f"coh = {ax_key[1]}")
            ax.set_ylim([0, 1.05])
            ax.set_ylabel('prop. right choices')
        elif 'PDW' in ax_key:
            ax.set_ylim([0, 1.05])
            ax.set_ylabel('prop. high bets')
        elif 'RT' in ax_key:
            # ax.set_ylim([0.5, 1.2])
            ax.set_ylabel('mean RT (s)')

        xhdgs = np.unique(data_obs['heading'])
        ax.set_xticks(xhdgs)
        ax.set_xticklabels(xhdgs, rotation=40, ha='right')

        if iax==0:
            ax.legend(handles=legend_handles, title=hue)

    # set overall xlabel at bottom of figure
    if hasattr(g.figure, 'supxlabel'):
        g.figure.supxlabel("Heading angle (°)")


    return g

# %%

def plot_rtq(
    RTq,
    row: Optional[str] = None,
    col: str = 'modality',
    hue: str = 'heading',
    depvar: str = 'PDW',
    palette = sns.color_palette(),
    **kwargs,
    ):

    g = sns.FacetGrid(
        RTq,
        row=row,
        col=col,
        hue=hue,
        aspect=1.5,
        height=6,
        sharey=True,
        palette=palette)

    # Iterate through each subplot and plot errorbars
    def plot_errorbars(**kwargs):
        plt.errorbar(
            x='RT_mean',
            y=depvar+'_mean',
            xerr='RT_cont_se',
            yerr=depvar+'_prop_se',
            **kwargs
            )
    g.map_dataframe(
        plot_errorbars, marker='.', capsize=3
        )

    # Customize labels and legend
    g.set_axis_labels('Mean RT (s)', f"Mean {depvar}")
    g.add_legend()

    return g


def RTquantiles(
    df: pd.DataFrame, 
    by_conds: list, 
    q_conds: Optional[list] = None, 
    nq: int=5, 
    depvar: str = 'PDW', 
    use_abs_hdg: bool = True
    ) -> pd.DataFrame:

    """
    Compute the quantiles of the RT and the dependent variable for each condition.
    Args:
        df (pd.DataFrame): the dataframe to compute the quantiles of
        by_conds (list): the conditions to group by
        q_conds (list): the conditions to group the quantiles by
        nq (int): the number of quantiles
        depvar (str): the dependent variable to compute the quantiles of
        use_abs_hdg (bool): whether to use the absolute heading
    Returns:
        pd.DataFrame: the dataframe with the quantiles of the RT and the dependent variable
    """

    q_conds = by_conds or q_conds

    if use_abs_hdg:
        df['heading'] = df['heading'].abs()

    # assign a quantile to each trial in the df, and store the mid of each quantile (for plotting)
    def calc_bin_edges(num_bins):
        q = np.arange(1, num_bins+1) / (num_bins+1)
        return np.concatenate(([0], q, [1]))

    transform_fcn = lambda x: pd.qcut(x, calc_bin_edges(nq), labels=False, duplicates='drop')
    df.loc[:, 'RTq'] = df.groupby(q_conds)['RT'].transform(transform_fcn)

    #qvals = df.groupby(q_conds)['RT'].transform(lambda x: pd.qcut(x, calc_bin_edges(x, nq)))
    #df.loc[:, 'qmid'] = qvals.apply(lambda x: x.mid)

    agg_funcs = {
        depvar: ['mean', prop_se, 'count'],
        'RT': ['mean', cont_se],
    }
    RTq = df.groupby(by_conds + ['RTq'])[[depvar, 'RT']].agg(agg_funcs).dropna(axis=0).reset_index()
    RTq.columns = ['_'.join(col) if col[0]==depvar or col[0]=='RT' else col[0] for col in RTq.columns]  # remove multi-level index

    return RTq

