%% dots3DMP human paper


% first, make sure we have cleaned datasets

% human non-RT
% human RT

% clean means no breakfix trials, only 2-target trials, only 'good'
% subjects, coherence and headings standardized, normalized conf
% see dots3DMP_cleanHumanData
% keeping raw cohs, headings and conf values as well

clear; clc; close all
cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/dataStructs
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes/'))


%%

subject = 'human';
conftask = 1;
RTtask   = 1; 


if ~RTtask
    load('human_20190625-20191231_nonRT_clean_Mar2022.mat');% human non-RT
    subjs = unique(data.subj);
    subjs2keep = 1:5;
%     subjs2keep = [1 3 4 5];
else
    load('human_20200213-20220317_RT_clean_Mar2022.mat') % human RT
    RTlims = [0.25 2.5];
    fnames = fieldnames(data);
    removethese = data.RT < RTlims(1) | data.RT > RTlims(2);
    for f=1:length(fnames), data.(fnames{f})(removethese) = []; end
    
    % should standardize RTs too for better analysis when pooling subject data?
    % let's check the distributions first
    
    subjs = unique(data.subj);
    subjs2keep = [1 3 5 8 9];
%     subjs2keep = [1 3 8 9];

end
fnames = fieldnames(data);
removethese = ~ismember(data.subj,subjs(subjs2keep));
for f=1:length(fnames), data.(fnames{f})(removethese) = []; end

% remove zero heading trials, insufficient data
% for f=1:length(fnames), data.(fnames{f})(data.heading==0) = []; end

% fix data.correct (see function description for issue!)
data.correct = dots3DMPCorrectTrials(data.choice,data.heading,data.delta);

% if doing pooled analysis, then just go forward with data
% otherwise re-assign data to be from a given subject and proceed the same
% way (except maybe we want to use raw conf values?)



%% basic data parsing and visualization

% these will be useful in multiple places
mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
% deltas = [-3 3]; % for showing only non-zero deltas
hdgs   = unique(data.heading);

% means per condition, logistic and gaussian fits
parsedData  = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit        = dots3DMP_fit_cgauss(data,mods,cohs,deltas,conftask,RTtask); 
% dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask) % eventually will just go straight to plotting fits of model?


%% 1. relationship between confidence and accuracy

% 1. conf vs P(right) for each modality
dots3DMP_plots_confchoice(parsedData,mods,cohs,deltas,conftask) 

%% 2. relationship between confidence and RT

% TODO review confidence distributions within each subject to determine
% threshold for high/lo (currently just the median normalized value)?

parsedData_byConf = dots3DMP_parseData_byConf(data,mods,cohs,deltas,hdgs,conftask,RTtask); 
gfit_byConf       = dots3DMP_fit_cgauss_byConf(data,mods,cohs,deltas,conftask,RTtask);
dots3DMP_plots_cgauss_byConf(gfit_byConf,parsedData_byConf,mods,cohs,deltas,hdgs,conftask,RTtask)

% 0 - plot errors/low bet only, 1 - plot correct/high bet only, 2 - plot correct/error or high/low bet separately, -1 - plot all trials
if RTtask
    dots3DMP_RTquantiles(data,conftask,-1); 
    confRT_distrs(data,mods,cohs,conftask,RTtask)
end


% Technical TODO
% make sure each figure is numbered and sized and saves, then make a
% wrapper that does analyses with each subject's data, by sub-sampling
% relevant trials from data



%% meta-sensitivity

% how good are subjects at reflecting on their own confidence
% across all
% dots3DMP_dprime_measures(data,mods,cohs,deltas) % NEED THE OUTPUTS!
%
% 
% 
% %% psychophysical cue weights
% % assume wvis always = 1-wves
% 
% wves = dots3DMP_cueWeights(gfit,cohs,deltas,conftask);
% 
% % bootstrapping for error bars, set rng for reproducibility
% rng(28);
% nboots = 20;
% [gfitBoot,wvesBoot] = dots3DMP_cgauss_bootstrap_func(data,gfit,mods,cohs,deltas,nboots,conftask,RTtask);
% dots3DMP_plotCueWeights(wves,wvesBoot,cohs,conftask) % plot the weights

