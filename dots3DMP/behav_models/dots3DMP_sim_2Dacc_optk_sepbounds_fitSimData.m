%% now try fitting the fake data to recover the generative parameters

clear
close all
load 2DAccSim_conftask2_162000trs

% convert to 0:1
if max(data.choice(~isnan(data.choice)))==2
    data.choice = logical(data.choice-1);    
end

% options.errfun = 'dots3DMP_fit_2Dacc_err_nSims';
options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noSim_CFtemp';
% options.nreps  = 100;
% options.confModel = 'evidence+time';

% SJ 10/2021, no longer doing model fits via Monte Carlo
options.runInterpFit = 0; 

options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% % % knoise = [0.07 0.07];
% % % sigmaVes = 0.03;
% % % sigmaVis = [0.03 0.03];

% initial guess (or hand-tuned params)
kmult = 65;
kvis  = kmult*cohs;
kves  = mean(kvis);
B = 1;
theta = 0.5;
Tnd     = [0.1 0.5 0.3]; % one for each modality
Ttc     = 0.35; % time to confidence!

if conftask==1 % SEP
    guess   = [kves kvis(1) kvis(2) B Tnd(1) Tnd(2) Tnd(3) Ttc];
    fixed   = [0    0       0       0 0      0      0      0  ];
else % PDW    
    guess   = [kves kvis(1) kvis(2) B Tnd(1) Tnd(2) Tnd(3) Ttc theta];
    fixed   = [0    0       0       0 0      0      0      0   0    ];
end

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


% check for fit-then-predict (missing values for +/- delta, etc)
parsedData = dots3DMP_parseData(fitInterp,mods,cohs,deltas,hdgs,conftask,RTtask); 

% plot it!
fitgauss = 1;
dots3DMP_plots_fit(data,fitInterp,conftask,RTtask,fitgauss);


