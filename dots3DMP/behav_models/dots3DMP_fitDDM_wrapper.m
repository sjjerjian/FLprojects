% dots3DMP_fitDDM_wrapper.m

% Generalized wrapper script for Dots DDM fitting, formerly a part of
% Dots_offlineAnalysis.m

% requires a struct data with at minimum a variable for choice and one for
% signed coherence

% CF updated 12/2021, again in 07/2022

clear; close all;


%% try fitting simulated data to recover the generative parameters

load tempsim.mat

%% or real data

load lucio_20220301-20221006_clean


%% some bookkeeping, then parse data, and plot if desired

if ~exist('allowNonHB','var'); allowNonHB=0; end

if allowNonHB==0
% ignore unabsorbed probability, ie no need for weighted sum with max_dur
% in RT calculation, e.g. when sim excluded trials that failed to hit bound
    options.ignoreUnabs = 1;
else
    options.ignoreUnabs = 0;    
end

options.RTtask = 1;
options.conftask = 2; % 1=continuous/rating, 2=PDW

% parse trial data into aggregated and other support vars
mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);
hdgs   = unique(data.heading);
RTCorrOnly = 0;
if ~exist('parsedData','var')  % e.g., if simulation was run
    parsedData = dots3DMP_parseData(data,mods,cohs,deltas,hdgs,options.conftask,options.RTtask);
end

% **** 
% optional [data will be plotted below regardless, along with the fits]
forTalk = 0;
% plot it
dots3DMP_plots(parsedData,mods,cohs,deltas,hdgs,options.conftask,options.RTtask)

% convert choice (back) to 0:1
if max(data.choice(~isnan(data.choice)))==2
    data.choice = logical(data.choice-1);    
end


%% now the fitting itself


%****** first select which model to fit ********
modelID=1; options.errfcn = @errfcn_DDM_2D_wConf_noMC; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (uses Wolpert's images_dtb_2d (method of images, from Moreno-Bote 2010))]
%***********************************************

% options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noMC';
% options.errfun = 'dots3DMP_fit_2Dacc_err_singlebound_noMC_unsigned';
options.errfun = 'dots3DMP_fit_2Dacc_err_singlebound_noMC_signed';

% options.nreps  = 100;
% options.confModel = 'evidence+time';

% SJ 10/2021, no longer doing model fits via Monte Carlo
options.runInterpFit = 1; 

options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% pass in orig params from sim as initial guess (or hand-tune if fixed=1)
if strfind(options.errfun,'singlebound') %%% paramNames = {'kves','kvisLo','kvisHi','B','TndVes','TndVis','TndComb','T2Conf','theta'};
    guess = [origParams.kves origParams.kvis(1) origParams.kvis(2) origParams.B/2 origParams.Tnds(1) origParams.Tnds(2) origParams.Tnds(3) origParams.ttc origParams.theta];
elseif strfind(options.errfun,'sepbounds') %%% paramNames = {'kves','kvisLo','kvisHi','BVes','BVis','BComb','muTnd','T2Conf','theta'};
    guess = [origParams.kves origParams.kvis(1) origParams.kvis(2) origParams.Bves origParams.Bvis origParams.Bcomb origParams.muTnd origParams.ttc origParams.theta];
end
% fixed = [0 0 0 0 0 0 0 0 0];
fixed = zeros(1,length(guess));

% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
fixed(:)=1;
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr  = 1;
options.RTtask   = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

if options.ploterr, options.fh = 400; end

%%
[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);
% fitInterp is in fact not obsolete, and needs fixing in ^^

% plot it!
dots3DMP_plots_fit_byCoh(data,fitInterp,conftask,RTtask);


%% in progress

% check for fit-then-predict (missing values for +/- delta, etc)

if any(unique(fit.delta(isnan(fit.choice)))==0) %--> should be nonzeros only
    keyboard
end
% this doesn't help, need to do the predict step in the fitDDM func




