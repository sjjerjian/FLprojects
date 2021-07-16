%% monte-carlo simulation of 1D drift-diffusion model w/ confidence

clear all; close all;

ntrials = 2000;

dbstop if error

cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';

dT = 1; % ms
maxdur = 2000;
timeAxis = 0:dT:maxdur;

k = 0.25; % drift rate coeff (conversion from %coh to units of DV)
B = 25; % bound height
mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % SD of momentary evidence
theta = 1; % threshold for Psure choice in units of log odds correct (Kiani & Shadlen 2009)

choice = nan(ntrials,1); % choices (left = -1, right = 1);
RT = nan(ntrials,1); % reaction time (or time-to-bound for fixed/variable duration task)
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);
pdw = nan(ntrials,1);

[logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,timeAxis,0);

%%
hitBound = zeros(ntrials,1);
logOddsCorr = zeros(ntrials,1);
tic
for t = 1:ntrials
    
    % the diffusion process
    dv = [0, cumsum(normrnd(mu(t),sigma,1,maxdur))];

    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(t) = maxdur;
        dvEnd = dv(RT(t));
        hitBound(t) = 0;
    else
        RT(t) = cRT;
        dvEnd = B*sign(dv(RT(t))); % cap bound crossings at bound value
        hitBound(t) = 1;
    end
   
    % choice
    if dvEnd==0
        choice(t) = sign(rand-0.5);
    else
        choice(t) = sign(dvEnd);
    end
    choice(t) = (choice(t)+1)/2; % convert to 0/1
    
    % look up the expected log odds corr for this trial's V and T
    whichV = find(abs(vAxis-dvEnd)==min(abs(vAxis-dvEnd)));
    whichT = find(abs(tAxis-RT(t))==min(abs(tAxis-RT(t))));   
    logOddsCorr(t) = logOddsCorrMap(whichV(1), whichT(1));

    % confidence rating
    expectedPctCorr(t) = logistic(logOddsCorr(t)); % convert to pct corr
    conf(t) = 2*expectedPctCorr(t) - 1; % convert to 0..1
    
    % post-decision wager:
    if logOddsCorr(t) >= theta % bet high
        pdw(t) = 1;
    elseif logOddsCorr(t) < theta % bet low
        pdw(t) = 0;
    else
        keyboard
    end    
    
end
toc

% add non-decision time (truncated normal dist)
TndMean = 300;
TndSD = 50;
TndMin = TndMean/2;
TndMax = TndMean+TndMin;
Tnd = zeros(ntrials,1);
for t = 1:ntrials
    while Tnd(t)<=TndMin || Tnd(t)>=TndMax % simple trick for truncating
        Tnd(t) = round(normrnd(TndMean,TndSD));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd;


% % quick sanity check that params are reasonable
% pCorrect_total = (sum(choice==1 & coh>0) + sum(choice==0 & coh<0)) / ntrials


%% plot proportion "rightward" (choice=1) and reaction time as a function of coherence

for c = 1:length(cohs)
    I = coh==cohs(c);
    pRight(c,1) = sum(I & choice==1) / sum(I);
    pHigh(c,1) = sum(I & pdw==1) / sum(I);
    meanRT(c,1) = mean(RT(I));
    meanConf(c,1) = mean(conf(I));
end

figure(121); clf; set(gcf,'Position',[900   465   900   750]);
subplot(2,2,1); plot(cohs,pRight(:,1),'co-'); ylim([0 1]);
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');

subplot(2,2,2); plot(cohs,pHigh(:,1),'bo-'); ylim([0 1]);
xlabel('Motion strength (%coh)'); ylabel('Proportion high bet');

subplot(2,2,3); plot(cohs,meanRT(:,1),'go-');
xlabel('Motion strength (%coh)'); ylabel('Reaction time (ms)');

subplot(2,2,4); plot(cohs,meanConf(:,1),'ro-'); ylim([0 1]);
xlabel('Motion strength (%coh)'); ylabel('Confidence rating');




