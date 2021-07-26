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
dots3DMP_plot_neuron_tuning(meanFRs(:,:,:,:,unit),semFRs(:,:,:,:,unit),cohs,hdgs);

%% Likelihood based decoding (Fig. 4)

step = 0.1;
numtrs = 100;
[tuning_curves, posterior, pop_lh, xq, simChoice] = ...
    dots3DMP_likelihood_decoding(data,meanFRs,hdgs,mods,cohs,deltas,numtrs,step);

%% Fig 4 in paper...

figure('position',[300 300 600 400],'color','w'); hold on
h = 5; % hdgs(5) is +1.2
d = 4;

% select a few simulated trials to plot
simTrs2plot = randperm(numtrs,5);

% vestibular only
ax=subplot(231); hold on;
                       %m,c,d,h, interphdgs, units
temp = squeeze(posterior(mods==1,cohs==16,deltas==0,h,:,simTrs2plot));
plot(xq,temp,'linew',1,'color','k')
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylab = ylabel('Likelihood (p(r|\theta))');
ylab.Position(2) = -0.1;
ylim([0 0.4]);

% visual only
ax=subplot(232); hold on;                    
plot(xq,squeeze(posterior(mods==2,cohs==60,deltas==0,h,:,simTrs2plot)),'color','r')
plot(xq,squeeze(posterior(mods==2,cohs==16,deltas==0,h,:,simTrs2plot)),'color','m','linestyle',':');
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylim([0 0.4]);
% combined, delta = +4, low coh
% visual is to the right of vestibular, at low coh, likelihood decoding is
% biased towards vestibular (left)
ax=subplot(234); hold on;                   
plot(xq,squeeze(posterior(mods==3,cohs==16,deltas==d,h,:,simTrs2plot)),'color','c','linestyle','-');
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylim([0 0.4]);
% combined, delta = +4, high coh
% visual is to the right of vestibular, at high coh, likelihood decoding is
% biased towards visual (right)
ax=subplot(235); hold on;                   
plot(xq,squeeze(posterior(mods==3,cohs==60,deltas==d,h,:,simTrs2plot)),'color','b','linestyle','-');
ax.XLim = [-6 6]; ax.TickDir = 'out'; ax.XTick = -5:2.5:5;
ylim([0 0.4]);
% #TODO - cumulative Gaussian fits to simChoices for psychometric curves
% (4e + f)
ax = subplot(233);
ax = subplot(236);
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
