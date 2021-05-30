% Kiani-style DDM w confidence

plotLOmap = 0; % plot log-odds corr map or not (set to 2 to see P(DV) also)


%% redefine pools if desired
rightpool = prefDirs<60 | prefDirs>300; % within 60 deg of zero
leftpool = prefDirs>120 & prefDirs<240; % within 60 deg of 180


%% calculate pooled spike counts and their difference

rightspikes = nan(nTrials,size(R,3));
leftspikes = nan(nTrials,size(R,3));

% Momentary evidence is the difference of spike counts
% (in a given window) between the two pools.
% First try every ms, although convolving with a 5-20 ms window is
% probably more realistic
windowWidth = 1;
convKernel = ones(1,windowWidth)/prod(windowWidth);
     
tic
% parfor_progress(nTrials);
    % parfor really struggles with this one for some reason
for t = 1:nTrials
    
    Rspikes = nanconv(sum(squeeze(R(t,rightpool,:))),convKernel,'same')'; % R = [trial, neuron, time]
    Lspikes = nanconv(sum(squeeze(R(t,leftpool,:))),convKernel,'same')';    
        
    % now the question is, over what time epoch to accumulate?
    % seems safe to ignore baseline; the trigger for the brain to start
    % accumulating is the dots onset, which we can assume is detected at
    % the same latency as MT responses.
%     Rspikes(1:latency) = []; % (note this shifts everything back in time by latency; this is
%     Lspikes(1:latency) = []; % fine because our models recover this time as non-decision time)
        % can save time by doing this after the loop is done
        
    % less clear is how much time to include after dots offset (e.g. the
    % rampdown) let's include only latency for now, such that the included
    % epoch is equal to duration
    Rspikes(dur(t)+latency+1:end) = NaN;
    Lspikes(dur(t)+latency+1:end) = NaN;
    
    rightspikes(t,:) = Rspikes;
    leftspikes(t,:) = Lspikes;
    
%     % take a look at the DV (accum. evidence) on a few trials:
%     if t<6
%         pooldiff = Rspikes-Lspikes;
%         figure(12); hold on; title(['coh = ' num2str(coh(t))]); plot(cumsum(pooldiff(t,:))); pause;
%     end

%     parfor_progress
end
toc

rightspikes(:,1:latency) = [];
leftspikes(:,1:latency) = [];


%% first simulate performance without PDW,
% to get params for the logOddsCorrect mapping

% % % dv = nan(nTrials,size(rightspikes,2)); % the DV
% % % dvAlt = dv; % temp, see below
choice = nan(nTrials,1); % vector to store choices
RT = nan(nTrials,1); % reaction time (or time-to-bound for fixed duration task)
finalV = nan(nTrials,1); % endpoint of the DV
hitBound = nan(nTrials,1);
Sigma = nan(nTrials,1);
Mu = nan(nTrials,1);

% need to guess bound and sigmaPooling; experiment with these to see their
% effect on choice+RT
B = 90;
sigmaPooling = 3;


tic
% parfor_progress(nTrials);
for t = 1:nTrials

    % add 'pooling noise', assuming the brain cannot pool and subtract
    % MT activity with arbitrary precision (Shadlen et al. 1996).
    % Practically speaking it's also necessary to get reasonable behavior
    % out of the simulation (poisson is way too good, and even fano=2 (high
    % end of realistic) is too good without pooling noise). Amount of noise
    % is somewhat arbitrary and depends on nNeurons and maxCorr.
    
    % unlike Shadlen et al. 1996, we're adding this at each time step, and 
    % to the difference variable, rather than to each pool's total spike
    % count separately (this should just affect the sigma needed)
    poolingNoise = round(sigmaPooling*randn(1,size(rightspikes,2)));
    
    diff = rightspikes(t,:) - leftspikes(t,:) + poolingNoise;
    dv = cumsum(diff);

    % save these for later
    Mu(t) = nanmean(diff);
    Sigma(t) = nanstd(diff);

    % temp: diagnosing mismatch between fittedk*C and actual Mu (need to
    % get a fittedk value (below) before trying this)
% % %     figure(t);plot(dv(t,:)); hold on;
% % %     dvAlt(t,1:dur(t)) = [0, cumsum(normrnd(fittedk*coh(t),Sigma(t),1,dur(t)-1))];
% % % %     dvAlt(t,1:dur(t)) = [0, cumsum(normrnd(Mu(t),Sigma(t),1,dur(t)-1))];
% % %     plot(dvAlt(t,:),'r-');
    
    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(t) = dur(t);
        finalV(t) = dv(RT(t));
        hitBound(t) = 0;
    else
        RT(t) = cRT;
        finalV(t) = B*sign(dv(RT(t))); % cap bound crossings at bound value
        hitBound(t) = 1;
    end

    if finalV(t)==0
        choice(t) = sign(rand-0.5);
    else
        choice(t) = sign(finalV(t));
    end
    
    choice(t) = (choice(t)+1)/2; % convert to 0/1
    
%     parfor_progress;
end
toc

% fit logistic
[beta, ~, ~] = glmfit(coh, choice, 'binomial');
xVals = min(coh):0.01:max(coh);
yVals = glmval(beta,xVals,'logit');

% from Shadlen 2006 book chapter:
% Pright =  1 / (1 + exp(-2kCB))
% beta = 2kB;
% thus,
fittedk = beta(2)/(2*B); % right?

% ^ this doesn't work. gives mu values >1 order of magnitude smaller than
% the actual (simulated) MT spikes would suggest. I leave this as an
% exercise for the reader to figure out why. :)


%% uncomment to plot & check for reasonable choice + RT

% add non-decision time (truncated normal dist)
TndMean = 300;
TndSD = 50;
TndMin = TndMean/2;
TndMax = TndMean+TndMin;
Tnd = zeros(nTrials,1);
for t = 1:nTrials
    while Tnd(t)<=TndMin || Tnd(t)>=TndMax % simple trick for truncating
        Tnd(t) = round(normrnd(TndMean,TndSD));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd;

for c = 1:length(cohs)
    I = coh==cohs(c);
    pRight(c,1) = sum(I & choice==1) /  sum(coh==cohs(c));
    RTmean(c,1) = mean(RT(I));
    RTse(c,1) = std(RT(I))/sqrt(sum(I));
    nc(c,1) = sum(I);
end
pRightSE = sqrt(pRight.*(1-pRight)./nc);

figure(2); clf; set(gcf,'Position',[900 320 480 760]);
subplot(2,1,1); 
errorbar(cohs,pRight,pRightSE,'o'); hold on;
plot(xVals,yVals); 
title('choice');
subplot(2,1,2); 
errorbar(cohs,RTmean,RTse,'s-');
title('RT');

% REMEMBER: in variable duration task, RT won't reflect realistic RT range,
% because:
proportion_hit_bound = sum(hitBound)/nTrials % is <1



%% now simulate PDW task

choice = nan(nTrials,1); % vector to store choices
pdw = nan(nTrials,1); % and wagers
RT = nan(nTrials,1); % reaction time (or time-to-bound for fixed duration task)
finalV = nan(nTrials,1); % endpoint of the DV
hitBound = nan(nTrials,1);
logOddsCorr = nan(nTrials,1);

% k = fittedk; % ignore (see above). Instead, use Mu calculated from spikes;
    % (near) zero c values are a problem for k=mu/c, need to use the dist
    % of Mu to map onto dist of k*c, or just omit the 0's?
I = abs(coh)> 0.01;
k = median(Mu(I)./coh(I));
    % median seems to work fine

% B = 30; % keep same as above
sigma = mean(Sigma);

% threshold (in log odds) for high bet (% hand tuned to get reasonable PDW)
theta = 1.5; 
% there's a weird mismatch between lookup table values of logOddsCorr and
% actual probability correct; see below (another exercise for the reader).
% in short, this theta is WAY too high.


% generate log odds corr map to implement PDW;
% set last argument to 0 to skip plotting (faster)
[logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,tAxis',plotLOmap);

tic
% parfor_progress(nTrials);
for t = 1:nTrials

    poolingNoise = round(sigmaPooling*randn(1,size(rightspikes,2)));
    
    diff = rightspikes(t,:) - leftspikes(t,:) + poolingNoise;    
    dv = cumsum(diff);
        
    % RT
    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(t) = dur(t);
        finalV(t) = dv(RT(t));
        hitBound(t) = 0;
    else
        RT(t) = cRT;
        finalV(t) = B*sign(dv(RT(t))); % cap bound crossings at bound value
        hitBound(t) = 1;
    end
   
    % choice
    if finalV(t)==0
        choice(t) = sign(rand-0.5);
    else
        choice(t) = sign(finalV(t));
    end
    
    % look up the expected log odds corr for this trial's V and T
    whichV = find(abs(vAxis-finalV(t))==min(abs(vAxis-finalV(t))));
    whichT = find(abs(tAxis-RT(t))==min(abs(tAxis-RT(t))));   
    logOddsCorr(t) = logOddsCorrMap(whichV(1), whichT(1));

    % post-decision wager:
    if logOddsCorr(t) >= theta % bet high
        pdw(t) = 1;
    elseif logOddsCorr(t) < theta % bet low
        pdw(t) = 0;
    else
        keyboard
    end
    
    choice(t) = (choice(t)+1)/2; % convert to 0/1
    
%    parfor_progress;
end
toc


% add non-decision time (truncated normal dist)
TndMean = 300;
TndSD = 50;
TndMin = TndMean/2;
TndMax = TndMean+TndMin;
Tnd = zeros(nTrials,1);
for t = 1:nTrials
    while Tnd(t)<=TndMin || Tnd(t)>=TndMax % simple trick for truncating
        Tnd(t) = round(normrnd(TndMean,TndSD));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd;



%% sanity check: choice probability (Britten et al. 1996)
I = abs(coh)<0.0001;
chooseright = I & choice==1; sum(chooseright);
chooseleft = I & choice==0; sum(chooseleft);
% first try based on full spike count
spCount = squeeze(nansum(R,3));
k=1;
for g = find(rightpool)
    CP(k,1) = rocN(spCount(chooseright,g), spCount(chooseleft,g), 50);
    k = k+1;
end
meanCP = mean(CP)







