% SJ 07-2021

% simple offline plots/analyses on behavioral data
% verify that new function codes all work


clear;clc

load('lucio_20210315-20210801_clean.mat')
RTtask = 1; 
conftask = 2;

% load('Human_20200213-20210707_clean.mat')
% RTtask = 1;
% conftask = 1;

% psychometric/chronometric plots

mods = unique(data.modality);
cohs = unique(data.coherence);
deltas = unique(data.delta); 
hdgs = unique(data.heading);

%% logistic fits

parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask);
dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask);
% dots3DMP_plots_func_forAumena(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask);

%% cum Gaussian fits... use gfit results for the weights and thresholds

gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask);
dots3DMP_plots_cgauss_func(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)


%% 
parsedDatabyConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask);
gfit = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask);


%% 

options.errfun = 'dots3DMP_fit_2Dacc_err_nSims';
options.nreps  = 100;
options.runInterpFit = 1; 

options.fitMethod = 'fms'; % 'fms','global','multi','pattern','bads'

options.confModel = 'evidence+time';
% choose whether to run fit with interpolated headings
% this is sort of redundant  for now, because model fits are
% generated via Monte Carlo and are going to be too noisy for a nice
% interpolated fit
options.runInterpFit = 1; 


options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% initial guess (or hand-tuned params)
kves    = 20;
kvis    = [15 35];
sigma   = [0.03 0.03 0.03];
BVes    = 1;
BVis    = 2;
BComb   = 1.5;
TndVes  = 300;
TndVis  = 300;
TndComb = 300;
fixed   = [0 1 1 1 1 1 1 1 1 1 1 1];
guess   = [kves kvis(1:2) sigma(1:3) BVes BVis BComb TndVes TndVis TndComb];

if conftask==2
    theta = 0.4;
    guess = [guess theta];
end
% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=1;
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr  = 1;
options.RTtask   = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

if options.ploterr, options.fh = 400; end

[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot data points and model fit curves
if options.runInterpFit 
    dots3DMP_plots_fit(data,fitInterp,conftask,RTtask)
end
