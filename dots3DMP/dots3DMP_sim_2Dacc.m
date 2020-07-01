% simulate "CCC" (cue-combination + confidence) experiment
% specifically in the case of the vis-ves heading task

% CF started it 2016
% heavily modified early 2019
% switched to 2D accumulator model 06/2020


% build expt and hand-pick some model params

clear
close all

plotExampleTrials = 0;

nreps = 500; % number of repetitions of each unique trial type
            % start small to verify it's working, then increase
            % (ntrials depends on num unique trial types)

cohs = [0.2 0.5]; % visual coherence levels (these are really just labels, since k's are set manually)
hdgs = [-10 -5 -2.5 -1.25 -eps eps 1.25 2.5 5 10]; % heading angles
deltas = [-3 0 3]; % conflict angle; positive means vis to the right
mods = [1 2 3]; % stimulus modalities: ves, vis, comb
duration = 2000; % stimulus duration (ms)

% sensitivity constants for converting heading angle into mean of momentary evidence
ks = 17; % scale factor, for quickly testing different k levels
  % set manually to get reasonable results from images_dtb_2d
kves = ks;
kvis = [.666*ks 1.5*ks]; % straddles vestibular reliability, by construction
    
sigmaVes = 0.03; % std of momentary evidence
sigmaVis = [0.03 0.03]; % allow for separate sigmas for condition, coherence
    % set manually to get reasonable looking dv trajectories;
    % also affects peakiness/flatness of confidence curve
    
B = 2; % bound height (also hand-tuned)
Tnd = 0.4; % non-decision time (eventually make this the mean of a dist)

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);
% % sigma = mean([sigmaVes sigmaVis]); % images_dtb_2d doesn't take a sigma argument

R.t = 0.001:0.001:duration/1000;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
R.lose_flag = 1;
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig

P =  images_dtb_2d(R);
% uses method of images to compute PDF of the 2D DV, as in van den Berg et
% al. 2016 (similar to Kiani et al. 2014)

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

% dv = cell(ntrials,1); % shouldn't need to store every trial's DV, but if you want to, it's here

choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1); % now this is the value of the losing accumulator
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

% for plotting
doneWith1 = 0;
doneWith2 = 0;
doneWith3 = 0;

% ME is now a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the flux file

tic
for n = 1:ntrials

    switch modality(n)
        case 1
            mu = kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
        case 2
            mu = kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
            % positive delta defined as ves to the left, vis to the right
            muVes = kves               * sind(hdg(n)-delta(n)/2) / 1000;
            muVis = kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;

            % optimal weights (Drugo et al.) 
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            mu = wVes.*muVes + wVis.*muVis;

            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb sigmaComb];
    end

    Mu = [mu,-mu]; % mean vector for 2D DV

    % convert correlation to covariance matrix
    V = diag(s)*S*diag(s);
    dv = [0 0; cumsum(mvnrnd(Mu,V,dur(n)-1))]; % bivariate normrnd
    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)

    % decision outcome
    cRT1 = find(dv(:,1)>=B, 1);
    cRT2 = find(dv(:,2)>=B, 1);
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000 + Tnd;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound, 
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000 + Tnd;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = dur(n)/1000 + Tnd;
            % this is interesting: which DV matters for confidence
            % if neither hits bound? for now take the (abs) maximum,
            % ie whichever was closest to hitting bound. (alternative
            % would be their average?)
        whichWon = dv(end,:)==max(dv(end,:));
        finalV(n) = dv(end,~whichWon); % the not-whichWon is the loser
        % % finalV(n) = mean(dvEnds);
        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000 + Tnd;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
        
    % use map to look up log-odds that the motion is rightward
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
    
    if plotExampleTrials
        if modality(n)==1 && hdg(n)==1.25 && choice(n)==1 && doneWith1==0 % make a better plot for talk/poster;
                              % must not shuffle trial list for this to go in correct in order
            figure(1000); set(gcf, 'Color', [1 1 1], 'Position', [100 100 350 375], 'PaperPositionMode', 'auto'); clf;
            
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'k-','LineWidth',2); hold on; 
            plot(1:length(dv),ones(1,length(dv))*B,'k-','LineWidth',4);
            xlabel('Time (ms)');
%             ylabel('Accum. evidence for rightward');
            ylim([-0.4*B B*1.1]); xlim([0 duration]);
%             set(gca,'yTick',-0.5:0.5:2);
            set(gca,'xTick',0:500:2000);
            changeAxesFontSize(gca,18,18);
            export_fig('ves_acc1','-eps');
            
            figure(1001); set(gcf, 'Color', [1 1 1], 'Position', [100 100 350 375], 'PaperPositionMode', 'auto'); clf;
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'r-','LineWidth',2); hold on; 
            plot(1:length(dv),ones(1,length(dv))*B,'k-','LineWidth',4);
            plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],[dv(round((RT(n)-Tnd)*1000),2) B],'r--')
            xlabel('Time (ms)');
%             ylabel('Accum. evidence for leftward');
            ylim([-0.4*B B*1.1]); xlim([0 duration]);
%             set(gca,'yTick',-0.5:0.5:2);
            set(gca,'xTick',0:500:2000);
            changeAxesFontSize(gca,18,18);
            export_fig('ves_acc2','-eps');

            doneWith1=1;
            n
        end
%         if modality(n)==2 && hdg(n)==2.5 && coh(n)==cohs(2) && choice(n)==1 && doneWith2==0
%             plot(dv,'r-','LineWidth',2); hold on; 
%             doneWith2 = 1;
%         end
%         if modality(n)==3 && hdg(n)==2.5 && coh(n)==cohs(2) && choice(n)==1 && delta(n)==0 && doneWith3==0
%             plot(dv,'b-','LineWidth',2); hold on;
%             doneWith3 = 1;
%             export_fig('threeTrials','-eps');
%         end
    end

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
data.RT = RT; % already in seconds
data.conf = conf;



%% plots

dots3DMP_parseData
dots3DMP_plots

% dots3DMP_parseData_splitConf
% dots3DMP_plots_splitConf



%% fit cumulative gaussians

dots3DMP_fit_cgauss


%% and plot them

dots3DMP_plots_cgauss



%% nicer looking versions

% dots3DMP_plots_cgauss_forTalk



% %% now try fitting the fake data to recover the generative parameters
%           %** AWAITS FITTING CODE FOR 2DACC
% 
% 
% % options.fitMethod = 'fms';
% % options.fitMethod = 'global';
% % options.fitMethod = 'multi';
% % options.fitMethod = 'pattern';
% options.fitMethod = 'bads';
% 
%     %    kves kvisMult B 
% fixed = [0    0        0];
% 
% % one small diff: in sim, kvis is just coh, here it will multiply coh
% 
% % initial guess (or hand-tuned params)
% kves = 1.2;
% kvisMult = 4; % will be multiplied by coh to get kvis (this simplifies parameterization)
% B = 70;
% 
% guess = [kves kvisMult B];
% 
% % ************************************
% % set all fixed to 1 for hand-tuning:
% % fixed(:)=1;
% % (can be used to fix some params and not others)
% % ************************************
% 
% % plot error trajectory (prob doesn't work with parallel fit methods)
% options.ploterr = 0;
% 
% [X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);
% 
% % plot it!
% %dots3DMP_plots_fit(data,fitInterp)
% 
% 



