%% simple analysis pipeline

clear;clc; close all
cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs

%%
subject = 'lucio';
conftask = 2; % 1=colorbars, 2=PDW
RTtask   = 1;
load('lucio_20210315-20210805_clean.mat') % recent lucio data, PDW + RT

% OR ...

% subject = 'human';
% conftask = 1; % 1=colorbars, 2=PDW
% RTtask   = 0;

% if ~RTtask, load('human_20190625-20191231_nonRT_clean.mat');% human non-RT, SfN 2021
% else,       load('human_20200213-20210922_RT_clean.mat') % human RT, SfN 2021
% end

%% basic parsing of data

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

% means per condition, logistic fits
parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 

%% gaussian fits
gfit = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 

%% plots
% separate subplots for each coh, with all mods on same
% or...separate subplots for each mod/delta, and all cohs on same

dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
% dots3DMP_plots_cgauss_byModDelta(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

%% psychophysical cue weights
% assume wvis always = 1-wves

wves = dots3DMP_cueWeights(gfit,cohs,deltas,conftask,1,0);

%% 
dots3DMP_RTquantiles(data,conftask,-1)

%%
nboots = 10;
[gfitBoot,wvesBoot] = ...
    dots3DMP_cgauss_bootstrap_func(data,gfit,mods,cohs,deltas,nboots);

%% confidence splits

parsedData_byConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit_byConf       = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask);

% plot it
dots3DMP_plots_cgauss_byConf(gfit_byConf,parsedData_byConf,mods,cohs,deltas,hdgs,conftask,RTtask)

%%
% weights separately for low and high coh?
% does this really make sense though - we don't know mapping between high conf in single
% cues to high conf in comb, or if one even exists
wves_byConf = dots3DMP_cueWeights(gfit_byConf,cohs,deltas,conftask,1);

% piecewise way
dots3DMP_ConfDelta(data,gfit,cohs,deltas,hdgs,conftask)
%% PDW and RT for correct vs incorrect trials
% specify mods and cohs to include only some, if desired

dots3DMP_CorrectVsErrorCurves(data,hdgs,mods,cohs,conftask,RTtask)
