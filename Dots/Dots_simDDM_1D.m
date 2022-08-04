% monte-carlo simulation of 1D drift-diffusion model w/ confidence
% formerly simDDM_simple; relies on Kiani 09 mapping
% CF started it circa 2016(?)


%% build expt and hand-pick some model params

clear all; close all;

ntrials = 50000;

dbstop if error

cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -eps eps 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';

dT = 1; % ms
max_dur = 2000;
tAxis = dT:dT:max_dur;

% params
k = 0.4; % drift rate coeff (conversion from %coh to units of DV)
B = 24; % bound height
mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % SD of momentary evidence
theta = 1.2; % threshold for high bet in units of log odds correct (Kiani & Shadlen 2009)

% % a very different set of params/outcomes, just to prove we can reproduce it in model code 
% k = 0.6;
% B = 40;
% mu = k*coh;
% sigma = 1;
% theta = 1;

theta2 = theta; % optional, a separate bet criterion for left v right
alpha = 0.4; % base rate of low bets (compresses & shifts the PDW curve down; does not interact w RT or choice)
TndMean = 300; % non-decision time (ms)
TndSD = 0; 
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



%% calculate log odds corr maps using Kiani 09 (FP4, Chang-Cooper) method

options.plot = 0; % plot marginal PDFs and LO map
makeLogOddsMap_1D


%% simulate drift-diffusion

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
    
    % the diffusion process
%     dv = [0, cumsum(normrnd(mu(n),sigma,1,max_dur))];
    dv = cumsum(normrnd(mu(n),sigma,1,max_dur));

    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(n) = max_dur;
        dvEnd = dv(RT(n));
        hitBound(n) = 0;
    else
        RT(n) = cRT;
        dvEnd = B*sign(dv(RT(n))); % cap bound crossings at bound value
        hitBound(n) = 1;
    end
   
    % choice
    if dvEnd==0 % if by some miniscule chance DV ends at zero, flip a coin
        choice(n) = sign(randn);
    else
        choice(n) = sign(dvEnd);
    end

    if hitBound(n)==0
        % look up the expected log odds corr for this trial's V and T
        thisV = find(abs(vAxis-dvEnd)==min(abs(vAxis-dvEnd)));
        thisT = find(abs(tAxis-RT(n))==min(abs(tAxis-RT(n))));
        if isnan(logOddsCorrMap(thisV(1), thisT(1))) % unfortunate rounding error: if dv gets very close to bound exactly at max_dur, but doesn't hit, it's effectively rounded up to the bound by thisV, but LO map is not defined at the bound
            if dvEnd>0
                logOddsCorr(n) = logOddsCorrMap(thisV(1)-1, thisT(1));
            else
                logOddsCorr(n) = logOddsCorrMap(thisV(1)+1, thisT(1));
            end
        else
            logOddsCorr(n) = logOddsCorrMap(thisV(1), thisT(1));
        end
        
        % post-decision wager:
        if logOddsCorr(n) >= theta % bet high
            pdw(n) = 1;
        elseif logOddsCorr(n) < theta % bet low
            pdw(n) = 0;
        else
            keyboard 
        end
        
    elseif hitBound(n)==1 % hit bound, no need for Pxt lookup
        if choice(n)==-1
            logOddsCorr(n) = log(Ptb_marginal(RT(n),1,1)./Ptb_marginal(RT(n),1,2));
            pdw(n) = bet_high_tb(RT(n),1);
        elseif choice(n)==1
            logOddsCorr(n) = log(Ptb_marginal(RT(n),2,2)./Ptb_marginal(RT(n),2,1));
            pdw(n) = bet_high_tb(RT(n),2);
        end        
    end
    
    % confidence rating
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
        
end
toc

simDDM_postProcessing % standardized for 1D and 2D

