% monte-carlo simulation of 1D drift-diffusion model w/ confidence
% formerly simDDM_simple; relies on Kiani 09 mapping
% CF started it circa 2016(?)


%% build expt and hand-pick some model params

clear all; close all;

ntrials = 5000;

dbstop if error

cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -eps eps 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';

dT = 1; % ms
maxdur = 2000;
timeAxis = 0:dT:maxdur;

% params
k = 0.25; % drift rate coeff (conversion from %coh to units of DV)
B = 25; % bound height
mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % SD of momentary evidence
theta = 1.1; % threshold for high bet in units of log odds correct (Kiani & Shadlen 2009)
TndMean = 300; % non-decision time
TndSD = 50; 
TndMin = TndMean/2;
TndMax = TndMean+TndMin;

origParams.k = k;
origParams.B = B;
origParams.sigma = sigma;
origParams.TndMean = TndMean;
origParams.TndSD  = TndSD;
origParams.TndMin = TndMin;
origParams.TndMax  = TndMax;



%% calculate log odds corr maps using Kiani 09 (FP4, Chang-Cooper) method

[logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,cohs,timeAxis,0);


%%

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
    dv = [0, cumsum(normrnd(mu(n),sigma,1,maxdur))];

    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(n) = maxdur;
        dvEnd = dv(RT(n));
        hitBound(n) = 0;
    else
        RT(n) = cRT;
        dvEnd = B*sign(dv(RT(n))); % cap bound crossings at bound value
        hitBound(n) = 1;
    end
   
    % choice
    if dvEnd==0
        choice(n) = sign(rand-0.5);
    else
        choice(n) = sign(dvEnd);
    end
    
    % look up the expected log odds corr for this trial's V and T
    whichV = find(abs(vAxis-dvEnd)==min(abs(vAxis-dvEnd)));
    whichT = find(abs(tAxis-RT(n))==min(abs(tAxis-RT(n))));   
    logOddsCorr(n) = logOddsCorrMap(whichV(1), whichT(1));

    % confidence rating
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
    
    % post-decision wager:
    if logOddsCorr(n) >= theta % bet high
        pdw(n) = 1;
    elseif logOddsCorr(n) < theta % bet low
        pdw(n) = 0;
    else
        keyboard
    end
    
end
toc

% add non-decision time (truncated normal dist)
Tnd = zeros(ntrials,1);
for n = 1:ntrials
    while Tnd(n)<=TndMin || Tnd(n)>=TndMax % simple trick for truncating
        Tnd(n) = round(normrnd(TndMean,TndSD));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd;

% quick sanity check that params are reasonable
pCorrect_total = sum(sign(choice)==sign(coh)) / ntrials

% should be >0.95 or else maxDur isn't long enough (or need urgency)
hitBoundPct = sum(hitBound)/length(hitBound)

%% format data as in experimental data files and generate output structs

coh(coh==0) = sign(randn)*eps; % should have no actual zeros, but if so, sign them randomly;
                               % this is just to assign a direction and correct/error
data.correct = choice==sign(coh);
data.direction = nan(ntrials,1);
data.direction(coh>0) = 0;
data.direction(coh<0) = 180;
coh(abs(coh)<1e-6) = 0; % now go back to one 'zero'
data.coherence = abs(coh);
data.scoh = coh;

data.choice = choice;
data.choice(data.choice==-1) = 0; % code elsewhere assumes 0s and 1s
data.RT = RT/1000; % convert to seconds
data.PDW = pdw;
data.conf = conf;

conftask = 2; % pdw (2) for now, even though conf rating can be generated here
RTtask = 1; RTCorrOnly = 0;
parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly);

% plot
cohs = unique(coh); wFit = 0; forTalk = 0;
Dots_plot(parsedData,cohs,conftask,RTtask,wFit,forTalk)



