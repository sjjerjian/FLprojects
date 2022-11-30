% Kiani-style DDM w confidence,
% operating on simulated (or real!) MT data, eg from simMT
% CF circa 2018-2020
clear;
% load("simMT_nNeu=360_nTr=1000.mat") % If running file on its own. As part of ms loop these vars will be provided.


%% Initialize PDW task variables


% generate log odds corr map to implement PDW;
% set last argument to 0 to skip plotting (faster)
if ms == 1 && tr == 1 % Only want to run these operations once
    addpath(genpath('FP4_Folder'));
    load('VarsFor_DDMDecode_PDW.mat')

    plotLOmap = 0; % plot log-odds corr map or not (set to 2 to see P(DV) also)
    choice = nan(nTrials,1); % vector to store choices
    pdw = nan(nTrials,1); % and wagers
    RT = nan(nTrials,1); % reaction time (or time-to-bound for fixed duration task)
    finalV = nan(nTrials,1); % endpoint of the DV
    hitBound = nan(nTrials,1);
    logOddsCorr = nan(nTrials,1);

    rightspikes = nan(nTrials,size(Ry,3));
    leftspikes = nan(nTrials,size(Ry,3));

    theta = 1.5; % threshold (in log odds), for high bet (% hand tuned to get reasonable PDW)
                 % there's a weird mismatch between lookup table values of logOddsCorr and
                 % actual probability correct; see below (another exercise for the reader).
                 % in short, this theta is WAY too high.

    I = abs(coh)> 0.01;
    k = median(Mu(I)./coh(I)); % "median seems to work fine"
    B = 90;       % Probly best to bring this in via VarsFor_DecodeDDM_PDW file.
    sigma = mean(Sigma);
    sigmaPooling = 3;

    [logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,tAxis',plotLOmap);

    rightspikes(:, 1:latency) = []; 
    leftspikes(:, 1:latency) = [];

elseif ms == 1 % Run once for each trial
    rightspikes(dur(t)+latency+1:end) = NaN; % Start vector after latency then fill times beyond RT with Nans
    leftspikes(dur(t)+latency+1:end) = NaN;
end

%% Run simulation

tic
% parfor_progress(nTrials);
for t = tr

    rightspikes(t,:) = sum(Ry(t,rightpool,1:ms),3,'omitnan'); % Fill in spike count totals for each pool, counts through till time ms
    leftspikes(t,:) = sum(Ry(t,leftpool,1:ms),3,'omitnan');

    poolingNoise = round(sigmaPooling*randn(1,size(rightspikes,2)));
    
    diff = rightspikes(t,1:ms) - leftspikes(t,1:ms) + poolingNoise;    
    dv = cumsum(diff); % cumulative difference between pool spCounts
% want to look up in logOddsMapL/R the value that corresponds to a
% particular dv (the dv closest to our current dv).
% Find index to the closest vAxis value to current dv.
% Enter this index and current ms as coordinates to locate prob value from
% logOddsMapR/L.

dv_lookup = find(min(abs(vAxis-dv(ms)))); % Find index to vAxis value closest to current dv, at time 'ms'
probL = logOddsMapL(dv_lookup,ms);
probR = logOddsMapR(dv_lookup,ms); % Will compare these values (maybe back in ms loop)
                                   % The bigger value indicates the more likely choice (RvL).
                                   % Identifying which neuron pool will be enhanced (RvL)
                                   % Can then apply modulation to L/Rpool
                                   % lambda values, as a function of probL/R


% logOddsMap plots dv values vs. time, so we can compare current 'dv' to
% find closest row value.
% What I want is the value associated with the color at the coordinates in
% plot of "Prob of DV for Rightward motion" and for leftward. Then can
% compare abs of these looked up values. The larger probability determines
% which pool will be enhanced. This is prob of R motion tho. NOT prob of R
% choice... how calc prob of R choice?

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
    whichV = find(abs(vAxis-finalV(t))==min(abs(vAxis-finalV(t)))); % Finds vAxis value closest to finalV. %YYY% seems redundant..? Why not just find(min(abs(vAxis-V(t))))
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
chooseright = I & choice==1; sum(chooseright); % total number of correct Rward trials
chooseleft = I & choice==0; sum(chooseleft); % total number of correct Lward trials
% first try based on full spike count
spCount = squeeze(sum(Ry,3,'omitnan'));
k=1;
for g = find(rightpool)
    CP(k,1) = rocN(spCount(chooseright,g), spCount(chooseleft,g), 50);
    k = k+1;
end
meanCP = mean(CP);


%% save the results
% Collect pdw, choice and npool variables to data file
% 
% disp('saving...');
% 
% clearvars -except hitBound dv finalV RT logOddsMapR logOddsMapL logOddsCorrMap vAxis pdw choice rightpool Ry tuning nNeurons nTrials prefDirs latency dir dur coh cohs poscohs tAxis
% % For feedback: logOddsMapR logOddsMapL logOddsCorrMap tAxis vAxis
% % For CP building: pdw choice rightpool Ry
% % For general: Ry tuning nNeurons nTrials prefDirs latency dir dur coh cohs poscohs tAxis
% tic
% save(sprintf('simMT_nNeu=%d_nTr=%d.mat', nNeurons, nTrials),'-v7.3');
% toc
% 
% disp('file saved.');







