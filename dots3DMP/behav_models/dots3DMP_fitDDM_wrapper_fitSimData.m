%% now try fitting the fake data to recover the generative parameters

clear
close all
load 2DAccSim_conftask2_90000trs_sigma05.mat
% load 2DAccSim_conftask2_180000trs_sigma04.mat

% convert to 0:1
if max(data.choice(~isnan(data.choice)))==2
    data.choice = logical(data.choice-1);    
end

% options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noMC';
% options.errfun = 'dots3DMP_fit_2Dacc_err_singlebound_noMC_unsigned';
options.errfun = 'dots3DMP_fit_2Dacc_err_singlebound_noMC_signed';

% options.nreps  = 100;
% options.confModel = 'evidence+time';

% SJ 10/2021, no longer doing model fits via Monte Carlo
options.runInterpFit = 0; 

options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% pass in orig params from sim as initial guess (or hand-tune if fixed=1)
if strfind(options.errfun,'singlebound') %%% paramNames = {'kves','kvisLo','kvisHi','B','TndVes','TndVis','TndComb','T2Conf','theta'};
    guess = [origParams.kves origParams.kvis(1) origParams.kvis(2) origParams.B origParams.Tnds(1) origParams.Tnds(2) origParams.Tnds(3) origParams.ttc origParams.theta];
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




