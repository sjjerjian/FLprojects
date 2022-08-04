% monte-carlo simulation of 2D drift-diffusion model (anticorrelated race)
% with confidence/PDW, for simulating simultaneous (peri-DW) or sequential
% (PDW) tasks

% CF spawned it from dots3DMP_sim_2Dacc_, merging with simDDM_1D_wConf
% (which was formerly simDDM_simple), in July 2022... surprised it didn't
% already exist (?)

%% build expt and hand-pick some model params

cd('/Users/chris/Documents/MATLAB/Projects/offlineTools/Dots')

clear all; close all
 
confModel = 'evidence+time';
% confModel = 'evidence_only';
% confModel = 'time_only';
% ^ these require different units for theta, time and evidence respectively

ntrials = 50000;

dbstop if error

cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -eps eps 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';

dT = 1; % ms
max_dur = 3000;

%% params
% these are very different from 1D model, k is larger and bound is smaller,
% maybe because of the way images_dtb is written (?)

k = 18; % drift rate coeff (conversion from %coh to units of DV)
B = 0.7; % bound height (0.6)
    mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % unit variance (Moreno-Bote 2010)           
theta = 1.2; % threshold for high bet in units of log odds correct (Kiani & Shadlen 09, 14)

% % alternate slate of params, to check robustness of fitting code (pre-param recovery)
% k = 8; % drift rate coeff (conversion from %coh to units of DV)
% B = 0.9;  % bound height (0.6)
%     mu = k*coh; % mean of momentary evidence (drift rate)
% sigma = 1; % unit variance (Moreno-Bote 2010)           
% theta = 1.6; % threshold for high bet in units of log odds correct (Kiani & Shadlen 09, 14)


alpha = 0; % base rate of low bets (offset to PDW curve, as seen in data)
TndMean = 300; % non-decision time (ms)
TndSD = 0; % 50-100 works well; set to 0 for fixed Tnd 
TndMin = TndMean/2;
TndMax = TndMean+TndMin;

origParams.k = k;
origParams.B = B;
origParams.sigma = sigma;
origParams.theta = theta;
origParams.alpha = alpha;
origParams.TndMean = TndMean;
origParams.TndSD  = TndSD;
origParams.TndMin = TndMin;
origParams.TndMax  = TndMax;

%%
% ME is a draw from bivariate normal with mean vector Mu and covariance
% matrix V start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the flux file

% convert correlation to covariance matrix
s = [sigma*sqrt(dT/1000) sigma*sqrt(dT/1000)]; % standard deviaton vector
         %^ variance is sigma^2*dt, so stdev is sigma*sqrt(dt), in seconds
V = diag(s)*S*diag(s); % equation for covariance matrix


%% calculate log odds corr map using Wolpert MoI method
tic
R.t = 0.001:0.001:max_dur/1000; % seconds
R.Bup = B;
R.drift = unique(abs(mu)); % unsigned drift rates
R.lose_flag = 1;
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig
P = images_dtb_2d(R);
toc

%% simulate bounded evidence accumulation

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
    ME = mvnrnd(Mu,V,max_dur-1);
    
% %     % check that it has the desired anticorr:
% %     CC = corrcoef(ME); corrcoef_me = CC(1,2)

    dv = [0 0; cumsum(ME)]; % bivariate normrnd
    
    % because Mu is signed according to direction (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision as in Kiani eqn. 3 and images_dtb

    % decision outcome
    cRT1 = find(dv(1:max_dur,1)>=B, 1); % time of first bound crossing(s):
    cRT2 = find(dv(1:max_dur,2)>=B, 1);
    
    % the possibilities are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound,
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = max_dur;
            % which DV matters for confidence if neither hits bound? 
            % SJ 07/2020 logOddsCorrMap is fixed, so just shift finalV up
            % so that 'winner' did hit bound,
        whichWon = dv(max_dur,:)==max(dv(max_dur,:));
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);
        % ^  shifting the losing dv up by whatever the
        % difference is between the bound and the winning dv

        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2]);
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
    
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
toc

simDDM_postProcessing % standardized for 1D and 2D




