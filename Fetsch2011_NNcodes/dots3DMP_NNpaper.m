% Using the dataset from Fetsch et al. 2011 Nature Neuroscience, to
% recreate figures/analyses

% SJ updated 06-2021

% some figures look slightly different from figures in paper...why?
% likelihood calculations (Fig. 4) seem to be off

%%

% load the data
cd /Users/stevenjerjian/Desktop/FetschLab/Analysis/data
load Fetsch_et_al_NatNeuro_2011.mat

% add relevant code folders to path, user specific
addpath(genpath('/Users/stevenjerjian/Desktop/FetschLab/Analysis/codes'))
addpath(genpath('/Users/stevenjerjian/Desktop/PhD/Codes/General/'))

%%
mods   = unique(data.modality);
cohs   = unique(data.coherence);
deltas = unique(data.delta);
hdgs   = unique(data.heading);

% select monkey (optional), [] = select both, 'W' = m18, 'Y' = m24
monkey = []; % [], 'W','Y'

newdata = data;

if ~isempty(monkey)
    switch monkey
        case 'W', monkID = m18;
        case 'Y', monkID = m24;
        otherwise
            error('no monkey with that ID');
    end
    removethese = ~startsWith(data.filename,['m' num2str(monkID)]);
    fnames = fieldnames(data);
    for F = 1:length(fnames)
        if strcmp(fnames(F), 'spikes')
            newdata.spikes(removethese,:) = [];
        else
            eval(['newdata.' fnames{F} '(removethese) = [];']);
        end
    end
end
%% Psychometric curves (Fig 1)

% logistic fits
parsedData = dots3DMP_parseData_func(newdata,mods,cohs,deltas,hdgs);
dots3DMP_plots_NN(parsedData,mods,cohs,deltas,hdgs);

% cum Gaussian fits... use gfit results for the weights and thresholds
gfit = dots3DMP_fit_cgauss_NN(newdata,mods,cohs,deltas);
dots3DMP_plots_cgauss_NN(gfit,parsedData,mods,cohs,deltas,hdgs)

%% weights and psycophysical thresholds (~Fig 2)
% NOTE: Fig. 2 in paper is plotted separately for each monkey

[wvesEmp,wvesPred] = dots3DMP_wgts_thres_NN(gfit.muPMF,gfit.sigmaPMF,cohs,deltas);

% bootstrapping for errorbars (resample data with replacement nboots times)
nboots = 100;
[muPMFboot,sigmaPMFboot,wvesEmpboot,wvesPredboot] = dots3DMP_cgauss_bootstrap_NN(newdata,mods,cohs,deltas,nboots);
dots3DMP_plot_wgts_bootstrap(wvesPred,wvesEmp,wvesEmpboot,wvesPredboot,sigmaPMFboot,gfit,cohs);

%% example MSTd neuron (Fig. 3)

[ufile,~,data.unitnum] = unique(data.filename,'stable');

clear monkUnit
monkUnit(startsWith(ufile,'m18'),1) = 'W'; % 48 units
monkUnit(startsWith(ufile,'m24'),1) = 'Y'; % 60 units

% calculate mean firing rate for each heading and condition
[meanFRs,semFRs] = dots3DMP_neuron_tuning(data,mods,cohs,deltas,hdgs); % all units

unit = 8; % 8 seems to be the one in the paper! 3 also looks nice
dots3DMP_plot_neuron_tuning(meanFRs,semFRs,unit,cohs,hdgs);

%% Likelihood based decoding (Fig. 4)

% define likelihood function, assuming independent Poisson variability
% tc = tuning curve of each neuron (interpolated)
% r_obs = observed population response
% lh = likelihood as function of heading

lh_func = @(tc,r_obs) exp(-tc) .* (tc .^ r_obs) ./ factorial(r_obs);

numunits = length(unique(data.unitnum));
numtrs   = 100; % number of simulated trials for each condition

% interpolate tuning curves at stepsize
step = 0.1;
xq = min(hdgs):step:max(hdgs);

% pre-allocate
tuning_curves = nan(length(mods),length(cohs),length(deltas),length(xq),numunits);
r_obs         = nan(length(mods),length(cohs),length(deltas),length(hdgs),numunits,numtrs);
lh            = nan(length(mods),length(cohs),length(deltas),length(hdgs),length(xq),numunits,numtrs);


% surely there's a way to vectorize this...anyway, loops ftw
for m=1:length(mods)
for c=1:length(cohs)
for d=1:length(deltas)
for u=1:numunits
    x = squeeze(meanFRs(m,c,d,:,u));
        
    if any(isnan(x)),continue,end % skip
 
    % tried using interp1, but causes errors... neurons (with low
    % firing rate) have the same exact firing rate for different headings?
    
    % linear regression
%     p = polyfit(hdgs,x,1);
%     y = polyval(p,xq);
%     y(y<0) = 0; % is this necessary?
    
    % or piecewise linear between headings % is this how Chris did it?  
    y = linear_pcw_fit(x,hdgs,0.1);

    tuning_curves(m,c,d,:,u) = y;
    
    % random draws of single-trial firing for each condition
    % ... and compute likelihood using tuning curve and function
    
    for h=1:length(hdgs)
        [~,idx] = min(abs(xq-hdgs(h)));
        r = poissrnd(tuning_curves(m,c,d,idx,u),numtrs,1);   

        r_obs(m,c,d,h,u,:) = r;
        
        for i=1:length(xq)
            tc = tuning_curves(m,c,2,i,u);
            lh(m,c,d,h,i,u,:) = lh_func(tc,r); 
        end
    end
end
end
end
end

%%

% SJ 06-2021 something is still off here (or in above calcs)...
% NOT YET WORKING

popn_llk = squeeze(prod(lh,6));         % product of llks across population
% can we use sum of logs here instead?

% normalize to obtain posterior (assuming uniform prior)
%popn_llk = popn_llk ./ max(sum(popn_llk,5),[],6);


% randomly select a few simulated trials
inds = randperm(size(popn_llk,6),10);

% plot example e.g. Fig 4a in paper
figure('position',[300 300 300 200]); ax = gca; hold on
% hdgs(5) is +1.2
                       %m,c,d,h, interphdgs, units
temp = squeeze(popn_llk(1,1,2,5,:,inds));

plot(xq,temp,'linew',1.5,'color','k')
ax.XLim = [-6 6];
ax.TickDir = 'out';
ax.XTick = -5:5:5;
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = -5:2.5:5;
% ax.YTick = 0:0.1:0.2;
% ax.YMinorTick = 'on';
% ax.YAxis.MinorTickValues = 0:0.05:0.2;
% ax.TickLength = [0.02 0.02];

%% Fig 4c,f --> simulated psychometric curves

% TO DO

%% Fig. 5 congruency index and tuning curves
% correlation between firing rate and heading for ves and vis separately,
% then product of these two

% use non-interpolated mean firing rates
vesFRs = squeeze(meanFRs(mods==1,1,2,:,:));
visFRs = permute(squeeze(meanFRs(mods==2,:,2,:,:)),[2 3 1]);
corrhdgs = hdgs;

% use interpolated tuning curves instead? no, probably artificially inflates
% correlations
% vesFRs = squeeze(tuning_curves(mods==1,1,2,:,:));
% visFRs = permute(squeeze(tuning_curves(mods==2,:,2,:,:)),[2 3 1]);
% corrhdgs = xq';

clear vesCorr visCorr* vesP visP
for u=1:size(vesFRs,2)
    [vesCorr(u,1),vesP(u,1)] = corr(vesFRs(:,u),corrhdgs);
    for c=1:length(cohs)
        [visCorr(u,c),visP(u,c)] = corr(visFRs(:,u,c),corrhdgs);
    end
end

[congruencyIndex,isSignificant] = dots3DMP_plot_congruency_NN(...
    monkUnit,tuning_curves,vesCorr,visCorr,vesP,visP,xq);

%% 
