temptemptemptemp



%% now try fitting the fake data to recover the generative parameters

if 0
% 
% % options.errfun = 'dots3DMP_fit_2Dacc_err_nSims';
options.errfun = 'dots3DMP_fit_2Dacc_err_sepbounds_noSim';

% % options.nreps  = 100;
% % options.confModel = 'evidence+time';


% choose whether to run fit with interpolated headings

% this is sort of redundant  for now, because model fits are
% generated via Monte Carlo and are going to be too noisy for a nice
% interpolated fit

% SJ 10/2021, no longer doing model fits via Monte Carlo
options.runInterpFit = 1; 

options.fitMethod = 'fms'; %'fms','global','multi','pattern','bads'
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% initial guess (or hand-tuned params)
kmult   = 30;
kvis    = kmult.*cohs';
kves    = mean(kvis);
BVes    = 0.6;
BVis    = 1.2;
BComb   = 0.8;
Tnd     = 300;
Ttc     = 0; % time to confidence!

fixed   = [0 1 1 1 1 1 1 1];
guess   = [kves kvis(1:2) BVes BVis BComb Tnd Ttc];

if conftask==2 % PDW
    theta = 0.6;

    fixed   = [0 1 1 1 1 1 1 1 1];

    guess   = [guess theta];
end

% ************************************
% set all fixed to 1 for hand-tuning, or 0 for full fit
fixed(:)=1;
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr  = 1;
options.dummyRun = 0;
options.RTtask   = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

if options.ploterr, options.fh = 400; end

% remove 
% removethese = data.RT == max(data.RT);
% fnames = fieldnames(data);
% 
% for f=1:length(fnames)
%     data.(fnames{f})(removethese) = [];
% end

[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot it!
dots3DMP_plots_fit_byCoh(data,fitInterp,conftask,RTtask,0)

end