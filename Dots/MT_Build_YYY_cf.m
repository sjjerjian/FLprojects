%% %YYY% 07/08/22 Temporally sequential MT_Build


%% build a simulated MT population and generate its spiking activity to a 
% standard set of random-dot motion stimuli
% CF 2014-2020

clear all
close all
    
%% build an experiment
    % Define number of trials to simulate and available trial stim coherences (cohs), Each trial is assigned a duration (dur) and stim coherence (coh) value

tic

% plot flags
plotMT = 1;
plotCorr = 1;
plotLOmap = 0; % for DDM plot log-odds corr map or not (set to 2 to see P(DV) also)


% noiseModel = 'indepPoisson';
% noiseModel = 'corrMat_Mazurek02';
noiseModel = 'corrMat_Shadlen96';

nTrials = 1200;

% signed motion coherence; negative is leftward
% cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -1e-8 1e-8 0.032 0.064 0.128 0.256 0.512];  % All possible coherence levels, signed
cohs = [-0.256 -0.064 -eps eps 0.064 0.256];  % All possible coherence levels, signed
coh = randsample(cohs,nTrials,'true')'; % Assign a coh value to each of N trials

dir = nan(size(coh));   % Create 'dir' vector 1 column, 1 row for each trial (nTrials = 1000)
dir(coh<0) = 180;   % Index to trials (row) with coh < 0. Assign these trials value 180 = Left is correct answer 
dir(coh>0) = 0;     % Same as above, but now Right is correct answer

% OR set dur to 1s (or any fixed value) -- btw makes it easier to verify fano factor (below)
dur(1:nTrials,1) = 1000;

toc

        disp('Task structure defined');
%% build MT population 
    % Each neuron assigned a preferred motion direction between 0,360 degrees
    % Figure2: Plot a neurons tuning curve response at all coh levels and all directions of motion
    % Figure3: Plots rightpool and leftpool frRates to al lstim directions at single coh level
    

tic

% nNeurons = 36;
% nNeurons = 180;
nNeurons = 360;

dirAxis = 0:(360/nNeurons):360-(360/nNeurons);
K = 3; % inverse variance term (will scale with coh)
ampl = 60;  % actual peak FR will be about half of this, at 51.2% coh
offset = 0; % this is not the same as spontaneous firing rate, which is of 
            % course >0; keeping this offset small or zero allows driven
            % firing rates to get close to zero for nonpreferred (null)
            % direction at high coherence, as MT does.
                % This has to do with baseline firing rates, to make sure
                % that fr rate goes down for anti-preferred directions of motion

prefDirs = linspace(dirAxis(1),dirAxis(end),nNeurons); % Each 'neuron' assigned a value between 0:360
poscohs = cohs(cohs>0); % Each unique pos coherence value
tuning = cell(1,length(poscohs)); % tuning holds 6 empty cells use {} to access
if plotMT   % Build figure outline
    figure(2); set(gca,'xlim',[0 360],'Xtick',0:60:360);
    xlabel('direction (deg)'); ylabel('firing rate (sp/s)');
    title('tuning curve for leftward preferring neuron (all cohs)');
    hold on;
end
for c = 1:length(poscohs)
    tuning{c} = zeros(nNeurons,length(dirAxis)); % Fill each cell in tuning with array of nNeurons x number of pref dirs (360 x 361) Each row in matrix corresponds to a particular neuron, each column a diff stim motion direction, the value occupying that spot is Fr rate response
    for n = 1:nNeurons
        % von Mises tuning functions (captures asymmetry in right/left response to coh)
        k = K * poscohs(c); % scale inverse variance by coh % k renewed value for each coherence
        tuning{c}(n,:) = ampl * exp(k*cosd(dirAxis-prefDirs(n))) / (2*pi*besseli(0,k)) + offset;  %YYYguess% Seems to calc firing rate of a single neuron across all stim motion directions (0:360)
                 %^ rows are neurons, columns are dirs           
    end
    if plotMT
        h(c) = plot(tuning{c}(90,:)); hold on; % Plot frRate values of neuron in row #180 (i.e. with prefDir 180) for all directions of motion % Store plot lines in h
        % h stores 6 lines, 1 for each coherence
%         Ltxt{c} = num2str(poscohs(c));
%         legend(h,Ltxt);
    end
end

% set baseline to average 0% coh response (or could be a bit lower);
% this will be used below during spike train generation
baseline = mean(mean(tuning{1}));

% Define two pools of neurons supplying evidence for the two alternatives;
% later can replace this with an arbitrary weight vector or smoothly
% varying weighting function (more realistic)
rightpool = prefDirs<60 | prefDirs>300; % within 60 deg of zero
leftpool = prefDirs>120 & prefDirs<240; % within 60 deg of 180


if plotMT
    figure(3); set(gcf,'Position', [300 300 1200 480],'Color',[1 1 1],'PaperPositionMode','auto'); clf; % Which of these does cool transparent plotlines?
    subplot(1,2,1); plot(dirAxis,tuning{end}(rightpool,:)); title('Right Pool'); % 0:360 (motion directions), Frate of all "rightpool" neurons, at cell 6 coh // stim direction v Neural response Fig 3. 
    xlabel('Motion direction (deg)'); ylabel('Firing rate (sp/s)');
    set(gca,'xtick',0:90:360,'xlim',[0 360],'tickdir','out','box','off');
    changeAxesFontSize(gca,20,20);
    subplot(1,2,2); plot(dirAxis,tuning{end}(leftpool,:)); title('Left Pool');
    xlabel('Motion direction (deg)'); ylabel('Firing rate (sp/s)');
    set(gca,'xtick',0:90:360,'xlim',[0 360],'tickdir','out','box','off');
    changeAxesFontSize(gca,20,20);
end
toc
        disp('MT neural population assembled');

%% assign spike counts based on tuning curve and noise model
    % Store all possible correlation values in vector 'c', turn into
    % a matrix 'Cor'with diag values manually set to 1

% create corr mat, if not independent poisson
if ~strcmp(noiseModel,'indepPoisson')
    maxCorr = 0.2;
    % pairwise correlation varies linearly from maxCorr to zero based on
    % tuning similarity (pref dir proximity)        %%YYY%% How can we
    % incorporate influence of task context on corr structure? i.e. max corr between neurons with similar pref Dir AND in same pool (right vs. left)

    fano = 1.8; % fano factor = variance/mean ratio (forced to be 1 for poisson,
                % but for real neurons it's often 1.5-2)
    c=linspace(maxCorr,0.0,nNeurons/2); c=[c fliplr(c)];    % c = CorrValue vector = 360 values mirrored at the 180:181 mid point .2:0,0:.2 % MORE Detail: distribute corr values across neural pop. do for one pool, then flip to make for second pool. Append second 180 corr values to first 180 values now have neurons 0:180 corr values = maxCorr:0, 181:360 = 0:maxCorr.
    Cor=toeplitz(c); % aka diagonal-constant matrix, shortcut for making the kind of corr mat we want
    Cor(logical(eye(size(Cor)))) = 1; % set main diagonal to be 1   %%YYY%% Logical converts the 1s and 0s of 'eye'dentity matric into indexible booleans :)
    tmpCor = Cor; tmpCor(tmpCor==1)=NaN; % for better color range when plotting
    if plotCorr
        figure(4); imagesc(tmpCor); set(gca,'ydir','normal'); axis square; colorbar;  % Corr with self = 1 along diagonal, 0 corr value band for neurons 180 degrees apart, symetric on both sides of central diag
        title('spike count correlation matrix: intended');
        xlabel('neuron ID (pref dir)'); ylabel('neuron ID (pref dir)');
    end
%     plotCorr = 0;
end

switch noiseModel

    % independent Poisson (simplest, but slow and unrealistic)
    case 'indepPoisson'
        Counts = nan(nTrials, nNeurons);
        tic
        for t = 1:nTrials
            c = poscohs == abs(coh(t));  % c = logical indexing to coh of trial t
            for n = 1:nNeurons
%                 lambda = RiD{c}(n,dir(t)==[0 180]) / 1000; % /1000 because tuning is in spikes/s and we're sampling every 1 ms
                lambda = tuning{c}(n,dirAxis==dir(t)) / 1000; % At coherence 'c', on trials with dir 't' (180 or 0), firing rate of neuron 'n' is /1000 because tuning response frRate is in spikes/s and we're sampling every 1 ms  %% tuning = matrix of coh specific frRates, n = a neuron with unique prefDir, dir = correct answer
                Counts(t,n) = sum(poissrnd(lambda, 1, dur(t)));  %% Lambda is defined by a neurons spProb per ms on trial 'tr'. Use Lambda to generate spCount vector for neuron n across all ms time bins of trial t. Sum across all counted spikes in trial t => trial spCount! store this value in Counts(t,n)
            end
        end
        toc

    % or, create spikes according to desired correlation matrix
    case 'corrMat_Mazurek02' % method of Mazurek et al. 2002 Nat Neurosci (intuitive, but slower)
        % % % for g = 1:20 % TEMP: to test that average corr value is as desired
        Counts = nan(nTrials, nNeurons);
        tic
        for t = 1:nTrials
            i = dir(t);
            c = poscohs == abs(coh(t));
%             m = RiD{c}(:,dir(t)==[0 180]); % mean
            m = tuning{c}(:,dirAxis==dir(t));  % frRate of all 360 neurons on trial with dir 't', at coh 'c'
            T = dur(t)/1000;
                % now we have a vector of mean counts and a correlation matrix based on
                % tuning similarity. So to get a covariance matrix just multiply by the
                % product of ith and jth s.d. where s.d. is sqrt(fano*m*T)
            sd = sqrt(fano*m*T);
            Cov = repmat(sd,1,nNeurons) .* repmat(sd,1,nNeurons)' .* Cor;
            Counts(t,:) = round(mvnrnd(m*T, Cov)); % multivariate normal random numbers
        end
        Counts(Counts<0)=0;
        toc
        % % % % plot pairwise responses for similar and unsimilar neurons, for a fixed stimulus (or not)
        % % % whichTr = abs(coh)<0.001;
        % % % durToNorm = dur(whichTr);
        % % % corrs1(g) = corr(Counts(whichTr,1)./durToNorm, Counts(whichTr,2)./durToNorm);
        % % % corrs2(g) = corr(Counts(whichTr,1)./durToNorm, Counts(whichTr,19)./durToNorm);
        % % % % figure; plot(Counts(whichTr,1),Counts(whichTr,2),'x'); title(num2str(corrs1(g)));
        % % % % figure; plot(Counts(whichTr,1),Counts(whichTr,19),'o'); title(num2str(corrs2(g)));
        % % % end
        % % % mean(corrs1)
        % % % mean(corrs2) % TEMP: to test that average corr value is as desired

    case 'corrMat_Shadlen96' % method of Shadlen et al. 1996 J Neurosci, Appendix 1 (>2x faster)
        tic
        rootQ = sqrtm(Cor);
        Counts = nan(nTrials, nNeurons);
        for t = 1:nTrials
            c = poscohs==abs(coh(t));
            T = dur(t)/1000;
%             m = RiD{c}(:,dir(t)==[0 180]);
            m = tuning{c}(:,dirAxis==dir(t));
            z = normrnd(0,1,nNeurons,1);
            y = rootQ * z;
            y = y.*sqrt(fano*m*T) + m*T; % variance is fano*mean, and SD is sqrt of that
            y(y<0) = 0; % no negative reponses
            Counts(t,:) = round(y); % Counts is a matrix of 't' trials by 'n' neurons, filled with the number of spikes counted on a trial, for each neuron
        end
        toc
        disp('Population correlation structure defined');
end

% below we'll also need spike *rates*, since dur varies
Rates = Counts./(dur/1000); % trials x neurons

%% check how well we achieved the intended corr mat & fano
% uncoment all if want to plot, then make plotFlag = 1
if plotCorr; figure(5); set(gcf,'position',[520 200 1500 1100]); end
Fanos = nan(nNeurons,length(cohs)); % cohs = pos and neg coh values % each row a neuron, each column a diff coh
for c = 1:length(cohs)
    I = coh==cohs(c); % Index to trials with coherence c
    ratemat = Rates(I,:);  % grab the spRate of all synthetic neurons on trials of coh 'c'
    corrmat = corrcov(cov(ratemat)); % Make covar matrix and turn into correlation matrix..

    %     % uncomment this to reassure yourself that corrcov() works:
    %     if c==1
    %     corrmat2 = nan(size(corrmat));
    %     for i = 1:nNeurons
    %         for j = 1:nNeurons
    %             corrmat2(i,j) = corr(ratemat(:,i),ratemat(:,j));
    %         end
    %     end
    %     figure;plot(corrmat(:),corrmat2(:),'x');
    %     end

    if plotCorr
        figure(5);
        subplot(3,4,c); imagesc(corrmat); set(gca,'ydir','normal'); axis square; colorbar;  % Plot correlation matrix
        xlabel('neuron ID (pref dir)'); ylabel('neuron ID (pref dir)');
        caxis([-0.1 maxCorr+0.1]); title(num2str(cohs(c)));
    end

    % variance to mean ratio across repeated presentations of the same
    % stimulus for all neurons; will be averaged across stimuli (cohs)
    % %% YYY Interpretation of above words... A neuron has its spRate variance calculated
    % across trials of repeated, identical, stimuli. The variance value
    % calculated for each unique stimulus are then averaged together to get
    % that neurons average variance.
    Fanos(:,c) = var(Counts(I,:))./mean(Counts(I,:));   % Index to all spCounts within a given coherence = Calc a coh wide FanoFactor... (not the interpretation reached from a"above words")
    %     Fanos(:,c) = var(Rates(I,:))./mean(Rates(I,:));   % Fanos has a value for each neuron, at each coh level
    % **this greatly overestimates the true Fano factor, even when using Rates,
    % because of added variance from variable durations. set dur(:)=1000 in
    % first code cell above to see the correct Fano

end
if plotCorr
    mean(Fanos) % Avg across all neurons
    mean(mean(Fanos))      %Avg across all cohs
    figure(5); title('spike count correlation matrix: actual (sep by coh)');
end


% confirm that variance of difference variable (momentary evidence) increases with coh, aka signal-dependent noise
%YYY% Interpretation... Difference in sp value across pools should grow as coh level grows.
%   Should variance in the difference value grow as well? Does poisson rate change for diff coh level? YES, it does! SO, bigger range of possible values of spRate, leads to more variance in spRate.

if plotMT
    rDiff = nan(nTrials,1);
    for t = 1:nTrials
        rDiff(t) = sum(Counts(t,rightpool)) - sum(Counts(t,leftpool));
    end
    for c = 1:length(poscohs)
        I = abs(coh)==poscohs(c); % Trials at coherence level 'c'
        rDiffMean(c) = mean(abs(rDiff(I)));
        rDiffVar(c) = var(rDiff(I));
    end
    figure(6); subplot(1,2,1); set(gcf,'position',[500 500 1000 420]);
    plot(poscohs,rDiffMean,'o-'); title('mean of difference variable (momentary evidence)'); xlabel('unsigned coh');
    subplot(1,2,2); plot(poscohs,rDiffVar,'v-'); title('variance of difference variable (momentary evidence)'); xlabel('unsigned coh');
end

%% Build SDF, Then generate spike vectors!
  
latency = 50; % ms
risetime = 20; % ms
dT = 1;
tAxis = 1:dT:max(dur)+3*latency; % it's 3x because 1x at the beginning, 2x at the end (ramp-down is slower than rise time)

SDF = nan(nTrials, nNeurons, length(tAxis));  % spProb for every neuron, on each trial, at every time point within that trial
Ry = nan(nTrials, nNeurons, length(tAxis));  % [trial x neuron x time(ms)] %Fill Ry with rastor vectors. Spike or no spike at each time point, of each trial, for each neuron. Leave as nans values past dur(tr) to avoid averaging into pop 

disp('Generating spike trains, may take tens of seconds to minutes...');

% first load the LO map, used to update spike rates during decision formation (ie feedback)
% must have been generated previously from vanilla simMT (and some assumptions)
load stimMT_LOmap_YYY.mat  % has some other vars in there too, including CP

% replot it if you like
% plotLOmap = 1;
% [logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,cohs,tAxis',plotLOmap);

sigmaPooling = 3; % forgot to add this to the .mat file

% we'll simulate behavior in this same loop, trial by trial, so initialize:
choice = nan(nTrials,1); % vector to store choices
pdw = nan(nTrials,1); % and wagers
RT = nan(nTrials,1); % reaction time (or time to covert bound for fixed duration task)
finalV = nan(nTrials,1); % endpoint of the DV
hitBound = nan(nTrials,1); % bound crossed or not
logOddsCorr = nan(nTrials,1); % log odds correct

tic  
for tr = 1:nTrials
    for n = 1:nNeurons
        meanrate = Rates(tr,n);   % Not sure why called "meanrate"... is just spRate for neuron 'n' on trial 'tr'
        sdf=[]; % spike density function
        sdf(1:latency) = baseline; % lay down some baseline rate, from tr=1 to vis latency
        sdf(end+1 : latency+risetime) = linspace(baseline,meanrate,risetime); % quick linear rise to mean
        decay = 1+0.2*sign(baseline-meanrate); % % steady state with some decay*, like real MT (up or down depending on whether mean rate is above or below baseline)
        sdf(end+1 : latency+dur(tr)) = linspace(meanrate, meanrate*decay, dur(tr)-risetime);
        sdf(end+1 : end+2*latency) = linspace(meanrate*decay, baseline, 2*latency); % include die-down over 2x latency

        % figure; plot(sdf); 

        % Add some NOISE to synthetic neural responses
        sdfNoiseSD = 0.2 * sdf; % arbitrary, but scales with spike rate
        sdfNoise = sdf + sdfNoiseSD.*randn(1,length(sdf));
        sdfNoise(sdfNoise<0)=0;

        % figure; plot(sdfNoise); 

        SDF(tr,n,1:length(sdf)) = sdfNoise; % Stores the spProb for each neuron, at each timepoint, on every trial.
    end

    trTime = dur(tr)+3*latency; % Time length of each trial
    DV = nan(1,trTime); % init DV

    poolingNoise = round(sigmaPooling*randn(1,trTime)); % generate this for each trial, but faster to do it outside the ms loop
    
    lambdaMod = zeros(1,nNeurons); % initialize to zero, will be updated each time through the loop based on DV (logOdds) in previous time step
    alpha = 0; % scaling factor for lambdaMod
    
    for ms = 1:dT:trTime
        lambda = SDF(tr,:,ms) / (1000/dT) + lambdaMod;
        Ry(tr,:,ms) = poissrnd(lambda); % Generate spike or no spike for time ms==1
        
        if ms==1
            DV(ms) = sum(Ry(tr,rightpool,ms)) - sum(Ry(tr,leftpool,ms)) + poolingNoise(ms);
        else
            DV(ms) = sum(Ry(tr,rightpool,ms)) - sum(Ry(tr,leftpool,ms)) + poolingNoise(ms) + DV(ms-dT);
        end
        
        % check for bound crossing
        if abs(DV(ms))>B
            hitBound(tr) = 1;
            choice(tr) = sign(DV(ms));
            RT(tr) = ms;
            finalV(tr) = B*sign(DV(ms)); % cap bound crossings at bound value

            % log odds correct is the value at the bound (not past it)
            DV_atB = sign(DV(ms))*B;
            whichV = find(abs(vAxis-DV_atB)==min(abs(vAxis-DV_atB)));
            whichT = find(abs(tAxis-ms)==min(abs(tAxis-ms)));            
            logOddsCorr(tr) = logOddsCorrMap(whichV(1), whichT(1));
            pdw(tr) = logOddsCorr(tr) >= theta;

            break; % ends the ms loop
             % (can't continue accumulating past bound because
             % LogOdds is not defined, thus can't update lambdaMod)
             
        elseif ms==trTime % out of time, bound not hit
            
        	hitBound(tr) = 0;
            choice(tr) = sign(DV(ms));
            RT(tr) = ms;
            whichV = find(abs(vAxis-DV(ms))==min(abs(vAxis-DV(ms))));
            whichT = find(abs(tAxis-ms)==min(abs(tAxis-ms)));            
            logOddsCorr(tr) = logOddsCorrMap(whichV(1), whichT(1));
            pdw(tr) = logOddsCorr(tr) >= theta;
            
        else % otherwise continue accumulating
            
            % map current DV to log odds
            whichV = find(abs(vAxis-DV(ms))==min(abs(vAxis-DV(ms))));
            whichT = find(abs(tAxis-ms)==min(abs(tAxis-ms)));   
            logOddsCorr_current = logOddsCorrMap(whichV(1), whichT(1));

            % for now, add spikes to the currently favored pool, but later
            % can try subtracting from the other one, or something else
            % entirely
            lambdaMod = zeros(1,nNeurons); % reinit
            if DV(ms)<0 % currently leftward is favored
                lambdaMod(leftpool) = 0.01 * alpha * logOddsCorr_current;
            elseif DV(ms)>0 % currently rightward is favored
                lambdaMod(rightpool) = 0.01 * alpha * logOddsCorr_current;
            end
        end                           % ^the 0.01 is a wild guess
        
    end
    
    % convince yourself DV looks reasonable
%     figure; plot(DV); pause

end
toc
disp('Spike vectors generated');

% add non-decision time (Tnd)
% TndMean = 300;
% TndSD = 50;
% TndMin = TndMean/2;
% TndMax = TndMean+TndMin;
% Tnd = zeros(nTrials,1);
% for tr = 1:nTrials
%     while Tnd(tr)<=TndMin || Tnd(tr)>=TndMax % simple trick for truncating
%         Tnd(tr) = round(normrnd(TndMean,TndSD));
%     end
% end

Tnd = 300;
DT = RT; % rename this 'decision time'
RT = DT+Tnd;

% convert choice to 0..1
if any(isnan(choice)); error; end
choice(choice==-1)=0;



%% Plot PDW Behavioural Data


simMT_plotResults


% For 'confirmation bias' fedback signal we want to selectively modify L
% and R pool firing rates (thus, future DV) in accordance with current
% state of DV. Log odds correct/L,R value measures state of DV.
%    If logOdds favors L, enhance lambda and vice versa for R.
% When we enhance the favored pool should the opposing pool be penalized?
% This may over amplify effect.. no? But is it experimentally supported?
% Dietrich 200.. microstim. % More likely to go right on weak right trials (and 0coh) AND less likely
% to go left on left motion trials (I think.. double check)






% %% Assesment plots of Neural spike trains
% 
%     % PSTHs, do they look reasonable?
% cohsForPSTH(:,1) = coh<-0.12; % strong left
% cohsForPSTH(:,2) = coh>-0.12 & coh<0.12; % weak
% cohsForPSTH(:,3) = coh>0.12; % strong right
% 
% clr = cool(size(cohsForPSTH,2));
% figure(8);
% minRT = 500;
% for c = 1:size(cohsForPSTH,2)
%     I = cohsForPSTH(:,c) & RT>=minRT;
%     sum(I)
%     plot(1:minRT-Tnd,mean(squeeze(mean(Ry(I,rightpool,1:minRT-Tnd),2,'omitnan')),'omitnan')*1000,'Color', clr(c,:)); hold on;
% %     plot(1:minRT,mean(squeeze(mean(Ry(I,rightpool,1:minRT),2,'omitnan')),'omitnan')*1000,'Color', clr(c,:)); hold on;
%                     % Ry = [trials, neurons, time]
% end
% title('right pool'); xlabel('time from dots onset (ms)'); ylabel('spikes/s');
% % legend(cellstr(num2str(cohsForPSTH')),'location','north','orientation','horizontal');
% 


%% save the results
% For use in simMT_test
% Make sure to initiate MT_test with nNeuron and nTrial values matching MT_build
% 
% disp('saving...');
% % 
% clearvars -except finalV hitBound logOddsCorrMap Mu Sigma pdw choice RT Tnd rightpool leftpool Ry nNeurons nTrials prefDirs latency dir dur coh cohs poscohs tAxis
% % 
% tic
% save(sprintf('simMT_nNeu=%d_nTr=%d.mat', nNeurons, nTrials),'-v7.3');
% toc
% % 
% disp('file saved.');



%% sanity check: choice probability (Britten et al. 1996)
% clear;
% load('simMT_nNeu=360_nTr=1000.mat')

I = abs(coh)<0.0001;
chooseright = I & choice==1;
sum(chooseright);
chooseleft = I & choice==0;
sum(chooseleft);

% first try based on full spike count
spCount = sum(Ry(:,:,latency+1:end),3,'omitnan');

% all cells, then select pools later
CP = nan(nNeurons,1);
for n = 1:nNeurons
    CP(n) = rocN(spCount(chooseright,n), spCount(chooseleft,n), 50);
end

meanCP = mean(CP)
meanCP_R = mean(CP(rightpool))
meanCP_L = mean(CP(leftpool))

figure;hist(CP(rightpool))
figure;hist(CP(leftpool))

%% sliding-window CPs

width = 50; % window width
step = 1; % window step (dt)
minRT = 700;
tAxisCP = width:step:minRT-Tnd;

CP = nan(nNeurons,length(tAxisCP));

for t = tAxisCP
    t
    spCount = sum(Ry(:,:,t-width+1:t),3,'omitnan');
    for n = 1:nNeurons
        CP(n,t) = rocN(spCount(chooseright & RT>=minRT,n), spCount(chooseleft & RT>=minRT,n), 50);
    end
end

%%
figure(2); clf;
set(gcf, 'Color', [1 1 1], 'Position', [200 500 750 450], 'PaperPositionMode', 'auto');
load CP00.mat
plot(tAxisCP,mean(CP(rightpool,tAxisCP)),'k-','LineWidth',3); hold on;
load CP02.mat
plot(tAxisCP,mean(CP(rightpool,tAxisCP)),'b-','LineWidth',3);
load CP04.mat; CP = CP+0.01;
plot(tAxisCP,mean(CP(rightpool,tAxisCP)),'r-','LineWidth',3);

set(gca,'ylim',[0.48 0.56],'ytick',0.48:0.02:0.56);
xlabel('Time (ms)');
ylabel('Choice probability');
changeAxesFontSize(gca,20,20);




