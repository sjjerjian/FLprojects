% monte-carlo simulation of 2D drift-diffusion model (anticorrelated race)
% with confidence/PDW, for simulating simultaneous (peri-DW) or sequential
% (PDW) tasks

% CF spawned it from dots3DMP_sim_2Dacc_, merging with simDDM_1D_wConf
% (which was formerly simDDM_simple), in July 2022... surprised it didn't
% already exist (?)

%% build expt and hand-pick some model params

clear; close all

datafolder = '/Users/chris/Documents/MATLAB';
codefolder = '/Users/chris/Documents/MATLAB/git/FLprojects/Dots/';

cd(codefolder)

confModel = 'evidence+time'; % default, Kiani 2014 etc
    % if want to test some alternatives:
% confModel = 'evidence_only'; 
% confModel = 'time_only';
    % ^ these require different units for theta, time and evidence respectively

ntrials = 50000;

dbstop if error

% set of coherences
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -eps eps 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';

% time step
dT = 1; % ms

% max duration, 
max_dur = 3000;

allowNonHB = 0; % allow non-hit-bound trials? if set to 0 and a trial lasts
% longer than max_dur, it is discarded. If set to 1, those trials are
% assigned RT = max_dur (affects comparison with mean RT in images_dtb, 
% which is calculated only for bound crossings)


%% PARAMS
% these are very different from 1D model: k is larger and bound is smaller,
% maybe because of the way images_dtb is written (?)

k = 18; % drift rate coeff (conversion from %coh to units of DV)
B = 0.7; % bound height
    mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % unit variance (Moreno-Bote 2010), not a free param           
theta = 1.2; % threshold for high bet in units of log odds correct (Kiani & Shadlen 09, 14)
alpha = 0.1; % base rate of low bets (offset to PDW curve, as seen in data)

% % alternate slate of params, to check robustness of fitting code (pre-param recovery)
% k = 8; % drift rate coeff (conversion from %coh to units of DV)
% B = 0.9;  % bound height
%     mu = k*coh; % mean of momentary evidence (drift rate)
% sigma = 1; % unit variance (Moreno-Bote 2010), not a free param             
% theta = 1.6; % threshold for high bet in units of log odds correct (Kiani & Shadlen 09, 14)
% alpha = 0.15; % base rate of low bets (offset to PDW curve, as seen in data)

 % Tnd = non-decision time (ms), to account for sensory/motor latencies
TndMean = 300;
TndSD = 0; % 50-100 works well; set to 0 for fixed Tnd 
TndMin = TndMean/2; % need to truncate the Tnd dist
TndMax = TndMean+TndMin;

% store the generative parameters, to use e.g. for (pre)param recovery
origParams.k = k;
origParams.B = B;
origParams.sigma = sigma;
origParams.theta = theta;
origParams.alpha = alpha;
origParams.TndMean = TndMean;
origParams.TndSD  = TndSD;
origParams.TndMin = TndMin;
origParams.TndMax  = TndMax;


%% calculate log odds corr map using Wolpert Method of Images code (van den Berg 2016, after Moreno-Bote 2010)
tic
R.t = 0.001:0.001:max_dur/1000; % seconds
R.Bup = B;
R.drift = unique(abs(mu)); % unsigned drift rates
R.lose_flag = 1;
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig
P = images_dtb_2d(R);
toc


%%
% momentary evidence is a draw from bivariate normal distribution
% with mean Mu (a two-vector) and covariance matrix V (2x2).

% Start with a desired correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the 'flux' file

% convert correlation to covariance matrix:
% s is the standard deviaton vector
s = [sigma*sqrt(dT/1000) sigma*sqrt(dT/1000)]; 
         %^ variance is sigma^2*dt, so stdev is sigma*sqrt(dt),
         % after converting from ms to seconds
         
V = diag(s)*S*diag(s); % equation for a covariance matrix


%% simulate bounded evidence accumulation

% preallocate
choice = nan(ntrials,1); % choices (left = -1, right = 1);
RT = nan(ntrials,1); % reaction time (or time-to-bound for fixed/variable duration task)
finalV = nan(ntrials,1); % now this is the value of the losing accumulator
hitBound = zeros(1,ntrials); % hit bound or not on that trial
logOddsCorr = nan(ntrials,1); % log odds correct
expectedPctCorr = nan(ntrials,1); % expected probability correct (converted to confidence rating)
conf = nan(ntrials,1); % confidence rating
pdw = nan(ntrials,1); % post-decision wager

tic
for n = 1:ntrials
    
    Mu = mu(n) * dT/1000; % mean of momentary evidence, scaled by delta-t
    Mu = [Mu; -Mu]'; % mean vector for 2D DV
    ME = mvnrnd(Mu,V,max_dur-1); % bivariate normrnd
    
% %     % uncomment to check that it has the desired anticorrelation:
% %     CC = corrcoef(ME); corrcoef_me = CC(1,2)

    dv = [0 0; cumsum(ME)];
    
    % Note: because Mu is signed according to direction (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision [as in Kiani eqn. 3 and images_dtb]

    % decision outcome
    cRT1 = find(dv(1:max_dur,1)>=B, 1); % time of first bound crossing(s):
    cRT2 = find(dv(1:max_dur,2)>=B, 1);
    
    % the possibilities are:
    % (1) only the 'right' accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only the 'left' accumulator hits bound,
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        hitBound(n) = 0;
        if allowNonHB
            RT(n) = max_dur;
                % which DV matters for confidence if neither hits bound? 
                % SJ 07/2020 logOddsCorrMap is fixed, so just shift finalV up
                % so that 'winner' did hit bound,
            whichWon = dv(max_dur,:)==max(dv(max_dur,:));
            finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);
            % ^  shifting the losing dv up by whatever the
            % difference is between the bound and the winning dv
            a = [1 -1];
            choice(n) = a(whichWon);
        else
            RT(n) = NaN;
            finalV(n) = NaN;
            choice(n) = NaN;
        end
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2]);
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end    
    
    if hitBound(n)==0 && allowNonHB==0
        logOddsCorr(n) = NaN;
        conf(n) = NaN;
        pdw(n) = NaN;
    else
        diffV = abs((P.y+B)-finalV(n));
        diffT = abs(R.t*1000-RT(n));
        switch confModel
            case 'evidence+time'
                % use map to look up log-odds that the motion is rightward
                thisV = find(diffV==min(diffV));
                thisT = find(diffT==min(diffT));
                logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
                expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
                conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
                pdw(n) = logOddsCorr(n) > theta;
            case 'evidence_only'
                conf(n) = max(diffV) ./ range(P.y);
                pdw(n) = max(diffV) > theta;
            case 'time_only'
                conf(n) = 1 - (RT(n)/1000) ./ range(P.t);
                pdw(n) = (RT(n)/1000) < theta;
        end
        if isnan(conf(n)), conf(n)=0; end % if dvs are almost overlapping, force conf to zero as it can sometimes come out as NaN
    end
    
end
toc

% run script for post-processing
simDDM_postProcessing % (standardized for 1D and 2D models)

