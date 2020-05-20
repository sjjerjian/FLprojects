% simulate "CCC" (cue-combination + confidence) experiment
% specifically in the case of the vis-ves heading task

% CF started it 2016
% heavily modified early 2019


% CCC emerged out of brainstorming my research plan for job apps and
% simultaneous discussions with John Morrison about perceptual confidence.
% The question is whether the estimation of sensory uncertainty that
% governs cue weighting also informs an explicit confidence report. Naively
% I imagined that such a result would support John's idea that perception 
% itself 'assigns' degrees of confidence (aka belief, aka subjective
% probability) -- the alternative being that confidence is solely post-
% perceptual. Resolving this particular debate, though of interest to
% philosophers of mind, is not really our primary goal although it can be
% useful (and fun) to engage with the arguments.

% The challenge for us is how to test whether weights and conf come from
% the same process, given that (a) they will share roughly the same mean
% even if not correlated from trial to trial (conditioned on the stimulus),
% and (b) confidence on multisensory trials isn't separable into components
% attributable to uncertainty of individual cues.

% Previously, I envisioned solving (b) by comparing empirical vs predicted
% weights, in two ways: (i) the traditional way based on single-cue
% thresholds, and (ii) a new way based on single-cue confidence ratings.
% This is similar to one of Devin's existing figures, but not exactly.
% Instead what I wanted to try is to calculate a predicted weight based on
% /relative/ confidence in the constituent single-cue conditions. How to
% calculate this? Presumably there's a normative definition. Let's find it.

% It should be noted that the way around this could be in the neural data.
% Decoding uncertainty from single-cue areas on single trials, leveraging
% what we learn from Miguel's experiment, should give a prediction for
% multisensory PDW...

% For now, let's just simulate a model with RT (Kiani 2014, or Drugo 2019)
% and see where it leads.

% ...

% Okay so the sim data look reasonable. What question can we ask with it?
% What's an alternative to test against? The question is, what is the
% nature of the weighting -- does it share signatures of the relationship
% between confidence and time, beyond what you'd expect from inverse
% variance as computed from single-cue choice data?
    % Analysis should be simple: compute weights for different RT quantiles
    % BUT, if both vis and ves confidence scales up/down equally w time,
    % their ratio won't change and hence neither will weights...
    
% How bout first just checking independent race model vs. integration?

%% build expt and hand-pick some model params

clear
close all

plotExampleTrials = 0;

model = 'kiani09+drugo14';
% model = 'indepRace';

nreps = 200; % number of repetitions of each unique trial type
            % start small to verify it's working, then increase
            % (ntrials depends on num unique trial types)

cohs = [0.2 0.5]; % visual coherence levels
hdgs = [-10 -5 -2.5 -1.25 -eps eps 1.25 2.5 5 10]; % heading angles
                            % (map fn seems to require even number of diff 
                            % levels (or just no zero), so we use +/- eps)
deltas = [-5 0 5]; % conflict angle; positive means vis to the right
mods = [1 2 3]; % stimulus modalities: ves, vis, comb
duration = 2000; % stimulus duration (ms)

% kves = 0.4; % sensitivity constant for converting heading angle into mean of momentary evidence
% kvis = [0.25 0.6]; % straddles vestibular reliability, by construction
% % vel and acc are now scaled 0..1, so scale up kves and kvis a bit to
% % compensate; kves and kvis would eventually be free parameters to be fit
% kves = kves*3;
% kvis = kvis*3;
%   moved this^ up here so it's easier to see/change the actual values (for param recovery)
kves = 1.2; % sensitivity constant for converting heading angle into mean of momentary evidence
kvis = [0.8 2]; % straddles vestibular reliability, by construction

B = 70; % bound height
sigmaVes = 1; % std of momentary evidence
sigmaVis = [1 1]; % allow for separate sigmas for condition, coherence
Tnd = 300; % non-decision time (or make this the mean of a dist?)

maxdur = duration;
% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);
sigma = mean([sigmaVes sigmaVis]);
[~, ~, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_3DMP(hdgs,k,B,sigma,maxdur,0);
% uses Fokker-Planck equation to propagate the probability density of the DV,
% as in Kiani & Shadlen 2009. Required for readout of confidence, although
% a simpler heuristic could be used (conf proportional to accum evidence)


% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020

% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
vel = normpdf(1:duration,duration/2,210);
vel = 0.37*vel./max(vel);
acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% figure; hh=plot(1:duration,vel,1:duration,acc);
% for i=1:length(hh)
%     hh(i).LineWidth=2;
% end
% ylim([-1.1 1.1]*max(acc))
% xlabel('time [s]')
% legend('Vel','Acc')
% changeAxesFontSize(gca,15,15);

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)

%% build trial list
% (can't just randsample the above vectors, because certain combinations of
% modality, coh, delta etc are invalid)
    
% repeat heading list once for ves, ncoh for vis, and ncoh x ndelta for comb
numHdgGroups = any(ismember(mods,1)) + ...
               any(ismember(mods,2)) * length(cohs) + ...
               any(ismember(mods,3)) * length(cohs)*length(deltas);
hdg = repmat(hdgs', numHdgGroups, 1);

lh = length(hdgs);
coh = nan(size(hdg));
delta = nan(size(hdg));
modality = nan(size(hdg));

% kluge for ves, call it the lowest coh by convention
if any(ismember(mods,1))
    coh(1:lh) = cohs(1); 
    delta(1:lh) = 0;
    modality(1:lh) = 1;
    last = lh;
else
    last = 0;    
end

if any(ismember(mods,2)) % loop over coh for vis
    for c = 1:length(cohs)
        these = last+(c-1)*lh+1 : last+(c-1)*lh+lh;
        coh(these) = cohs(c);
        delta(these) = 0;
        modality(these) = 2;
    end
    last = these(end);
end

if any(ismember(mods,3)) % loop over coh and delta for comb
    for c = 1:length(cohs)
        for d = 1:length(deltas)
            here = last + 1 + (c-1)*lh*length(deltas) + (d-1)*lh;
            these = here:here+lh-1;
            coh(these) = cohs(c);
            delta(these) = deltas(d);
            modality(these) = 3;
        end
    end
end

% now replicate times nreps and shuffle (or not):
condlist = [hdg modality coh delta];
% trialTable = repmat(condlist,nreps,1); 
trialTable = Shuffle(repmat(condlist,nreps,1),2);
    % why shuffle? well, sometimes it's easier to find particular trial
    % types to test out when debugging
hdg = trialTable(:,1);  
modality = trialTable(:,2);  
coh = trialTable(:,3);  
delta = trialTable(:,4);  
ntrials = length(trialTable);


% sample durations from truncated exponential?
    % not practical experimentally.
    % instead, make constant dur, assume RT task, or (equivalently, as far as
    % the simulation is concerned) fixed duration with internal bounded
    % process (Kiani et al. 2008)
dur = ones(ntrials,1) * duration;

%% bounded evidence accumulation

% assume momentary evidence is proportional to sin(heading),
% as in drugowitsch et al 2014

% dv = nan(ntrials,max(dur)); % shouldn't need to store every trial's DV, but if you want to, it's here
% if strcmp(model,'indepRace')
%     dvVes = dv;
%     dvVis = dv;
% end

choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1);
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

% for plotting
doneWith1 = 0;
doneWith2 = 0;
doneWith3 = 0;


tic
for n = 1:ntrials

%     % slow version: codebank(71)
    
    % faster version: avoids FOR loop over the variable t
    switch modality(n)
        case 1
            mu = acc .* kves * sind(hdg(n)); % mean of momentary evidence
%             dv = [0, cumsum(normrnd(mu,sigmaVes,1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaVes))];
        case 2
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n));
%             dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n)),1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n))))];

        case 3
            % positive delta defined as ves to the left, vis to the right
            muVes = acc .* kves               * sind(hdg(n) - delta(n)/2);
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n) + delta(n)/2);
            
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );

            mu = wVes.*muVes + wVis.*muVis;
            % here the DV is a sample from a dist with mean = weighted sum
            % of means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
%             dv = [0, cumsum(normrnd(mu,sigmaComb,1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaComb))]; % if mu is a vector

            if strcmp(model,'indepRace')
%                 dvVes = [0, cumsum(normrnd(muVes,sigmaVes,1,dur(n)-1))];
%                 dvVis = [0, cumsum(normrnd(muVis,sigmaVis(cohs==coh(n)),1,dur(n)-1))];
                
                dvVes = [0, cumsum(normrnd(muVes,sigmaVes))];
                dvVis = [0, cumsum(normrnd(muVis,sigmaVis(cohs==coh(n))))];
            
            end                
    end
        

    if plotExampleTrials
        if modality(n)==1 && hdg(n)==2.5 && doneWith1==0 % make a better plot for talk/poster;
                              % must not shuffle trial list for this to go in correct in order
            figure(1000); set(gcf, 'Color', [1 1 1], 'Position', [100 100 450 350], 'PaperPositionMode', 'auto'); clf;
            plot(dv,'k-','LineWidth',2); hold on; 
            plot(1:length(dv),ones(1,length(dv))*B,'k-','LineWidth',4);
            plot(1:length(dv),ones(1,length(dv))*-B,'k-','LineWidth',4);
            xlabel('Time (ms)'); ylabel('Accum. evidence (DV)');
            ylim([-B-1 B+1]); xlim([0 2000]);
            set(gca,'yTick',-70:35:70,'Yticklabel',{'-70' '-35' '0' '35' '70'});
            % (evidence is shown continuing to accumulate
            % past the bound, although it realy stops there)
            xlabel('Time (ms)'); ylabel('Accumulated evidence');
            changeAxesFontSize(gca,20,20);
            doneWith1=1;
        end
        if modality(n)==2 && hdg(n)==2.5 && coh(n)==cohs(2) && doneWith2==0
            plot(dv,'r-','LineWidth',2); hold on; 
            doneWith2 = 1;
        end
        if modality(n)==3 && hdg(n)==2.5 && coh(n)==cohs(2) && delta(n)==0 && doneWith3==0
            plot(dv,'b-','LineWidth',2); hold on;
            doneWith3 = 1;
%             export_fig('threeTrials','-eps');
        end
    end
    
    
    if strcmp(model,'indepRace') && modality(n)==3
        cRTves = find(abs(dvVes)>=B, 1);
        cRTvis = find(abs(dvVis)>=B, 1);
        % the options are:
        % (1) only ves hits bound,
        if ~isempty(cRTves) && isempty(cRTvis)
            RT(n) = cRTves + Tnd;
            finalV(n) = B*sign(dvVes(cRTves));
            hitBound(n) = 1;
        % (2) only vis hits bound, 
        elseif isempty(cRTves) && ~isempty(cRTvis)
            RT(n) = cRTvis + Tnd;
            finalV(n) = B*sign(dvVis(cRTvis));
            hitBound(n) = 1;
        % (3) neither hits bound,
        elseif isempty(cRTves) && isempty(cRTvis)
            RT(n) = dur(n) + Tnd;
                % this is interesting: which DV matters for confidence
                % if neither hits bound? for now take the (abs) maximum,
                % ie whichever was closest to hitting bound. (alternative
                % would be their average?)
            dvEnds = [dvVes(end) dvVis(end)];
            whichWon = abs(dvEnds) == max(abs(dvEnds));
            finalV(n) = dvEnds(whichWon);
%             finalV(n) = mean(dvEnds);
            hitBound(n) = 0;
        % (4) or both do
        else
            RT(n) = min([cRTves cRTvis]) + Tnd;
            dvEnds = [dvVes(cRTves) dvVis(cRTvis)];
            whichWon = [cRTves<=cRTvis cRTves>cRTvis];
            finalV(n) = B*sign(dvEnds(whichWon));
            hitBound(n) = 1;
        end
    else
        cRT = find(abs(dv)>=B, 1);
        if isempty(cRT) % did not hit bound
            RT(n) = dur(n) + Tnd;
            finalV(n) = dv(dur(n));
            hitBound(n) = 0;
        else % hit bound
            RT(n) = cRT + Tnd;
            finalV(n) = B*sign(dv(cRT));
            hitBound(n) = 1;
        end    
    end
    choice(n) = sign(finalV(n));
    
    % use map to look up log-odds that the motion is rightward
    diffV = abs(vAxis-finalV(n));
    diffT = abs(tAxis-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = logOddsCorrMap(thisV(1), thisT(1));
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1

end
toc

choice(choice==0) = sign(randn); % not needed under usual circumstances

% sanity check:
pCorrect_total = (sum(choice==1 & hdg>0) + sum(choice==-1 & hdg<0)) / ntrials



%% format data like the real expt

choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

conftask = 1;
data.modality = modality;
data.heading = hdg;
data.coherence = coh;
data.delta = delta;
data.choice = choice;
data.RT = RT;
data.conf = conf;



%% plots

dots3DMP_parseData
dots3DMP_plots

dots3DMP_parseData_splitConf
dots3DMP_plots_splitConf

%% fit cumulative gaussians
% (needed for weights calculation)

cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
    % for probabilities, error is negative log likelihood of observing the data, which is
    % [ log(Pright(hdg)) + log(1-(~Pright(hdg))) ]
cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice))))); 

flippedGauss = @(b,hdg) 1 - ( min(max(b(1),0),1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));
    % for continuous values, error is sum squared error
flippedGauss_err = @(param,SEP,hdg) sum((flippedGauss(param,hdg)-SEP).^2);

% unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)
% dots3DMP_plots_cgauss


%%



%% now try fitting the fake data to recover the generative parameters


% options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
options.fitMethod = 'bads';

    %    kves kvisMult B 
fixed = [0    0        0];

% one small diff: in sim, kvis is just coh, here it will multiply coh

% initial guess (or hand-tuned params)
kves = 1.2;
kvisMult = 4; % will be multiplied by coh to get kvis (this simplifies parameterization)
B = 70;

guess = [kves kvisMult B];

% ************************************
% set all fixed to 1 for hand-tuning:
% fixed(:)=1;
% (can be used to fix some params and not others)
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr = 0;

[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot it!
%dots3DMP_plots_fit(data,fitInterp)





