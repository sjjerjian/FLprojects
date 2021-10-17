%% dots3DMP behavioral plots and analyses
% SJ 10/2021

% for SfN Poster

% datasets
% lucio - RT, Mar to Aug 2021
% human - non-RT up to Feb 2020
% human - RT, Mar 2020 to Oct 2021
% simulated data, different models

clear; clc; close all
cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs

%% select subject, load the data

subject = 'human';

switch subject
    
    case 'lucio'
        load('lucio_20210315-20210805_clean.mat') % recent lucio data, PDW + RT
        conftask = 2; % 1=colorbars, 2=PDW
        RTtask   = 1;
        
    case 'human'

        conftask = 1;
        RTtask   = 1; % change this to select RT or non-RT dataset
       
        if ~RTtask, load('human_20190625-20191231_nonRT_clean.mat');% human non-RT, SfN 2021
        else,       load('human_20200213-20210922_RT_clean.mat') % human RT, SfN 2021
        end
        
    case 'simul' % load simulated data

        load('2DAccSim_conftask2_8100trs.mat')
        % conftask & RTtask should already be saved in file
end


%% basic parsing of data

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

% means per condition, logistic fits
parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 

% gaussian fits
gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 

%% plots

% logistic fits
% dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

% separate subplots for each coh, with all mods on same subplot
dots3DMP_plots_cgauss_byCoh(gfit,parsedData,1,cohs,deltas,hdgs,conftask,RTtask)

% or...separate subplots for each mod/delta, and all cohs on same subplot
% this one needs tidying to look nice
% dots3DMP_plots_cgauss_byModDelta(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

%% psychophysical cue weights
% assume wvis always = 1-wves

wves = dots3DMP_cueWeights(gfit,cohs,deltas,conftask);

% bootstrapping for error bars
nboots = 100;
[gfitBoot,wvesBoot] = dots3DMP_cgauss_bootstrap_func(data,[],mods,cohs,deltas,nboots,0,RTtask);

% plot the weights
dots3DMP_plotCueWeights(wves,wvesBoot,cohs,0)

%% Confidence as function of decision time (RT quantiles)

if ~isfield(data,'correct')
    data.correct = (data.heading>0 & data.choice==2) | (data.heading<0 & data.choice==1) | (data.heading==0 & rand<0.5);
end

% third argument specifies which trials to use 
% -1: all, 0: error trials only, 1: correct trials only
% -1 will also lead to plotting of p(correct) as function of RT quantiles
% 0 will only use weaker stimuli (set to bottom 3), since stronger stimuli
% produce few errors 

dots3DMP_RTquantiles(data,conftask,-1)

%% PDW and RT for correct vs incorrect trials

dots3DMP_CorrectVsErrorCurves(data,conftask,RTtask,1)

%% Psychometric curves split by PDW

% Choice curves for low bet should be shallower - indicative that PDW is
% predictive of accuracy, even within a stimulus condition!

parsedData_byConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit_byConf       = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask);

% plot it
dots3DMP_plots_cgauss_byConf(gfit_byConf,parsedData_byConf,mods,cohs,deltas,hdgs,conftask,RTtask)

%% relationship between confidence and cue weights

% calculate weights separately for low and high coh?
% CF notes - this probably doesn't make sense to try and do however - 
% we don't know mapping between high conf in single cues to high conf in comb, or if one even exists!

% anyway, it's easy to do technically, just use gfit_byConf from above
    
% wves_byConf = dots3DMP_cueWeights(gfit_byConf,cohs,deltas,conftask,1);

% or we can do a piecewise way (part 2 of this function)
% this function actually runs three separate analyses relating confidence
% and cue conflict

% 1. compare shifts of choice mu vs conf mu in non-conflict and conflict conditions 
% 2. piecewise comparison of PRight for fixed absolute heading, different conflicts as function of coh
% 3. average confidence in conflict vs no conflict (low headings only)

dots3DMP_ConfDelta(data,gfit,cohs,deltas,hdgs,conftask)


%% MODELLING

options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noSim';
options.runInterpFit = 1; 

options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'

% initial guess (or hand-tuned params)
kves    = 25;
kmult   = 50;
kvis    = kmult.*cohs';
BVes    = 0.9;
BVis    = 1.5;
BComb   = 1.1;
Tnd     = 300;
Ttc     = 300; % time to confidence!

guess   = [kves kvis(1:2) BVes BVis BComb Tnd Ttc];

if conftask==2 % PDW
    theta = 0.8;
    guess = [guess theta];
end

% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
fixed(:)=0;
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr  = 1;
options.dummyRun = 0;
options.RTtask   = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

if options.ploterr, options.fh = 400; end

[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot it!
dots3DMP_plots_fit_byCoh(data,fitInterp,conftask,RTtask,0);

