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

subject = 'lucio';

switch subject
    
    case 'lucio'
        load('lucio_20210315-20210805_clean.mat') % recent lucio data, PDW + RT
        conftask = 2; % 1=colorbars, 2=PDW
        RTtask   = 1;
        
        RTlims = [0.25 2.25];

    case 'human'

        conftask = 1;
        RTtask   = 1; % change this to select RT or non-RT dataset
       
        if ~RTtask
            load('human_20190625-20191231_nonRT_clean.mat');% human non-RT, SfN 2021
            
        else
            load('human_20200213-20210922_RT_clean.mat') % human RT, SfN 2021
            RTlims = [0.25 2.5];
        end
        
        
    case 'simul' % load simulated data

        load('2DAccSim_conftask2_8100trs.mat')
        % conftask & RTtask should already be saved in file
end

if RTtask
    fnames = fieldnames(data);
    removethese = data.RT < RTlims(1) | data.RT > RTlims(2);
    
    for f=1:length(fnames)
        data.(fnames{f})(removethese) = [];
    end
end

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

if strcmp(subject,'human') && RTtask
    hdgs(hdgs==0) = []; % not enough good data, so let's just remove for now
end

%% basic parsing and summary plots of data

% means per condition, logistic fits
parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 

% gaussian fits
gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 

% logistic fit plots
% dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

% separate subplots for each coh, with all mods on same subplot
dots3DMP_plots_cgauss_byCoh(gfit,parsedData,1,cohs,deltas,hdgs,conftask,RTtask)

% or...separate subplots for each mod/delta, and all cohs on same subplot
% n.b. - this one needs tidying to look nice
% dots3DMP_plots_cgauss_byModDelta(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

%% psychophysical cue weights
% assume wvis always = 1-wves

wves = dots3DMP_cueWeights(gfit,cohs,deltas,conftask);

% bootstrapping for error bars
rng(28); % for reproducibility
nboots = 5;
[gfitBoot,wvesBoot] = dots3DMP_cgauss_bootstrap_func(data,gfit,mods,cohs,deltas,nboots,conftask,RTtask);

% plot the weights
dots3DMP_plotCueWeights(wves,wvesBoot,cohs,conftask,gfit)

%% Confidence as function of decision time (RT quantiles)

if ~isfield(data,'correct')
    data.correct = (data.heading>0 & data.choice==2) | (data.heading<0 & data.choice==1) | (data.heading==0 & rand<0.5);
end

% third argument specifies which trials to use 
% -1: all, 0: error trials only, 1: correct trials only
% -1 will also lead to plotting of p(correct) as function of RT quantiles
% 0 will only use weaker stimuli (set to bottom 3), since stronger stimuli
% produce few errors 


% TODO new version which simplifies plotting of correct or error, high or
% low bet
% dots3DMP_RTquantiles(data,conftask,-1)
dots3DMP_RTquantiles2(data,conftask,1)

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

dots3DMP_ConfDelta(data,gfit,cohs,deltas,hdgs,conftask,RTtask,3)


%% MODELLING

% options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noSim';
options.errfun = 'dots3DMP_fit_2Dacc_err_noSim';
options.fitMethod = 'multi'; %'fms','fmsbnd','global','multi','pattern','bads'
options.whichFit = 3; % 0 - choice only, 1 - choice + RT, 2 - choice + conf, 3 - ALL
options.paramNames = {'kves','kvisLo','kvisHi','BVes','BVis','BComb','TndVe','TndVi','TndCo','T-Conf','theta','cLapse'};

fixAll = 0;

% initial guess (or hand-tuned params)
% kmult   = 0.3;
% kvis    = kmult.*cohs';
% kves    = mean(kvis);
% kves    = 46;
% B       = 0.5;
% Tnd     = [0.4 0.5 0.4];
% Ttc     = 0; % time to confidence, ignored for now...

kves = 0.2;
kvis = [0.1 0.3];
B    = [0.3 0.3 0.3];
Tnd  = [0.5 0.65 0.55];
Ttc  = 0;
cL   = [0.1]; % confLapse rate, lapse rate of high bets

guess   = [kves kvis B Tnd Ttc];

fixed = [0 0 0 0 0 0 0 0 0 1];


if conftask==2 % PDW
    theta = 0.08;
    guess = [guess theta cL];
    fixed = [fixed 0 0];
end

if options.whichFit == 0
    fixed = [0 0 0 0 0 0 1 1 1 1 1 1];
end

options.lowerbound = guess/2;
options.upperbound = guess*2;


% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
if fixAll, fixed(:)=1; end
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr  = 1;
options.dummyRun = 0;
options.RTtask   = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

if options.ploterr, options.fh = 400; end

options.runInterpFit = 1;
[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);
dots3DMP_plots_fit_byCoh(data,fitInterp,conftask,RTtask,0);
    
