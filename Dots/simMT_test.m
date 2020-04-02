% given MT population responses (real or simulated), simulate decision/pdw
% tasks under different candidate mechanisms for constructing confidence
% CF 2014-2020

clear; close all;

% simulated MT data will have been saved with filename defined by nNeurons
% and nTrials
nNeurons = 360;
nTrials = 1000;

disp('loading...');
tic
load(sprintf('simMT_nNeu=%d_nTr=%d.mat', nNeurons, nTrials));
toc

ustim = zeros(size(coh));

%% temp/optional: simulate microstimulation

stimpool = prefDirs<=10 | prefDirs>=350; % within 10 deg of zero
lambda = 0.02; % think of this as spikes added per ms [doesn't take much!]
tic
for t = 1:nTrials
    if round(rand)==1 % random half of trials
        win = latency+1:latency+dur(t);
        Rtemp = squeeze(R(t,stimpool,win));
        Rtemp = Rtemp + poissrnd(lambda, 1, dur(t));
        R(t,stimpool,win) = Rtemp;
        ustim(t) = 1;
    end
end
toc


% now choose an option:
%% DDM (Kiani & Shadlen, uses FP4.mex)
disp('running DDM...');
simMT_decode_DDM
disp('done');
simMT_plotResults

if sum(ustim)>0; simMT_plotResults_ustim; end


%% classic PPC or likelihood re/decoding (Ma et al., Jazayeri & Movshon)
disp('running PPC...');
simMT_decode_PPC
disp('done');

% % temp: quickly try different gammas
% gamma = 10;
% pdw(betaAll(:,2)>=gamma) = 1;
% pdw(betaAll(:,2)<gamma) = 0;
% simMT_plotResults

if sum(ustim)<10
    simMT_plotResults
else
    simMT_plotResults_ustim;
end
% no shift in confidence curve for ustim, because it's only based on
% variance of posterior! obviously this doesn't work...


%% sampling (Hoyer & Hyvarinen, Fiser et al. (Haefner, Berkes))
disp('running sampling...');
simMT_decode_sampling
disp('done');

% % temp: quickly try different gammas
% gamma = 0.4;
% pdw(betaAll(:,2)>=gamma) = 1;
% pdw(betaAll(:,2)<gamma) = 0;
% simMT_plotResults

if sum(ustim)<10
    simMT_plotResults
else
    simMT_plotResults_ustim;
end


%% gain variability (Goris)

% this requires some modification to simMT_build...




