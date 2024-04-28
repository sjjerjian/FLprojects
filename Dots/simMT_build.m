% build a simulated MT population and generate its spiking activity to a 
% standard set of random-dot motion stimuli
% CF 2014-2020

clear
close all
    
%% build an experiment

tic

% plot flags
plotMT = 0;
plotCorr = 0;

% noiseModel = 'indepPoisson';
% noiseModel = 'corrMat_Mazurek02';
noiseModel = 'corrMat_Shadlen96';

nTrials = 1000;

% signed motion coherence; negative is leftward
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -1e-8 1e-8 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,nTrials,'true')';

dir = nan(size(coh));
dir(coh<0) = 180;
dir(coh>0) = 0;

% viewing duration drawn from truncated exponential distribution
dur = zeros(nTrials,1);
meandur = 360; mindur = 80; maxdur = 900;
for n = 1:nTrials
    while dur(n)<=mindur || dur(n)>=maxdur % simple trick for truncating
        dur(n) = round(exprnd(meandur));
    end
end

% OR set dur to 1s (or any fixed value) -- btw makes it easier to verify fano factor (below)
dur(1:nTrials,1) = 1000;

toc

%% build MT population
tic

% nNeurons = 36;
% nNeurons = 180;
nNeurons = 360;

dirAxis = 0:360/nNeurons:360-360/nNeurons;

K = 3; % inverse variance term (will scale with coh)
ampl = 60;  % actual peak FR will be about half of this, at 51.2% coh
offset = 0; % this is not the same as spontaneous firing rate, which is of 
            % course >0; keeping this offset small or zero allows driven
            % firing rates to get close to zero for nonpreferred (null)
            % direction at high coherence, as MT does

prefDirs = linspace(dirAxis(1),dirAxis(end),nNeurons);
poscohs = cohs(cohs>0);
tuning = cell(1,length(poscohs));
if plotMT
    figure(2); set(gca,'xlim',[0 360],'Xtick',0:60:360);
    xlabel('direction (deg)'); ylabel('firing rate (sp/s)');
    title('tuning curve for leftward preferring neuron (all cohs)');
    hold on;
end
for c = 1:length(poscohs)
    tuning{c} = zeros(nNeurons,length(dirAxis));
    for n = 1:nNeurons
        % von Mises tuning functions (captures asymmetry in right/left response to coh)
        k = K * poscohs(c); % scale inverse variance by coh
        tuning{c}(n,:) = ampl * exp(k*cosd(dirAxis-prefDirs(n))) / (2*pi*besseli(0,k)) + offset;
                 %^ rows are neurons, columns are dirs           
    end
    if plotMT
        h(c) = plot(tuning{c}(round(nNeurons/2),:)); hold on;
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
    figure(3); set(gcf,'Position', [300 300 1200 480],'Color',[1 1 1],'PaperPositionMode','auto'); clf;
    subplot(1,2,1); plot(dirAxis,tuning{end}(rightpool,:)); title('Right Pool');
    xlabel('Motion direction (deg)'); ylabel('Firing rate (sp/s)');
    set(gca,'xtick',0:90:360,'xlim',[0 360],'tickdir','out','box','off');
    changeAxesFontSize(gca,20,20);
    subplot(1,2,2); plot(dirAxis,tuning{end}(leftpool,:)); title('Left Pool');
    xlabel('Motion direction (deg)'); ylabel('Firing rate (sp/s)');
    set(gca,'xtick',0:90:360,'xlim',[0 360],'tickdir','out','box','off');
    changeAxesFontSize(gca,20,20);
end

toc

%% assign spike counts based on tuning curve and noise model

% create corr mat, if not independent poisson
if ~strcmp(noiseModel,'indepPoisson')
    maxCorr = 0.2;
    % pairwise correlation varies linearly from maxCorr to zero based on
    % tuning similarity (pref dir proximity)
    fano = 1.8; % fano factor = variance/mean ratio (forced to be 1 for poisson,
                % but for real neurons it's often 1.5-2)
    c=linspace(maxCorr,0.0,nNeurons/2); c=[c fliplr(c)];
    Cor=toeplitz(c); % aka diagonal-constant matrix, shortcut for making the kind of corr mat we want
    Cor(logical(eye(size(Cor)))) = 1; % set main diagonal to be 1
    tmpCor = Cor; tmpCor(tmpCor==1)=NaN; % for better color range when plotting
    if plotCorr
        figure(4); imagesc(tmpCor); set(gca,'ydir','normal'); axis square; colorbar;
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
            c = poscohs == abs(coh(t));
            for n = 1:nNeurons
%                 lambda = RiD{c}(n,dir(t)==[0 180]) / 1000; % /1000 because tuning is in spikes/s and we're sampling every 1 ms
                lambda = tuning{c}(n,dirAxis==dir(t)) / 1000; % /1000 because tuning is in spikes/s and we're sampling every 1 ms
                Counts(t,n) = sum(poissrnd(lambda, 1, dur(t)));
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
            m = tuning{c}(:,dirAxis==dir(t));
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
            Counts(t,:) = round(y);
        end
        toc

end

% below we'll also need spike *rates*, since dur varies
Rates = Counts./(dur/1000);

%% check how well we achieved the intended corr mat & fano

if plotCorr; figure(5); set(gcf,'position',[520 200 1500 1100]); end
Fanos = nan(nNeurons,length(cohs));
for c = 1:length(cohs)
    I = coh==cohs(c);
    ratemat = Rates(I,:);
    corrmat = corrcov(cov(ratemat));
    
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
        subplot(3,4,c); imagesc(corrmat); set(gca,'ydir','normal'); axis square; colorbar;
        xlabel('neuron ID (pref dir)'); ylabel('neuron ID (pref dir)');
        caxis([-0.1 maxCorr+0.1]); title(num2str(cohs(c)));
    end
    
    % variance to mean ratio across repeated presentations of the same
    % stimulus for all neurons; will be averaged across stimuli (cohs)
    Fanos(:,c) = var(Counts(I,:))./mean(Counts(I,:));
%     Fanos(:,c) = var(Rates(I,:))./mean(Rates(I,:));
% **this greatly overestimates the true Fano factor, even when using Rates, 
% because of added variance from variable durations. set dur(:)=1000 in 
% first code cell above to see the correct Fano

end
if plotCorr
    mean(Fanos)
    mean(mean(Fanos))
    figure(5); suptitle('spike count correlation matrix: actual (sep by coh)');
end


% confirm that variance of difference variable (momentary evidence)
% increases with coh, aka signal-dependent noise
if plotMT
rDiff = nan(nTrials,1);
for t = 1:nTrials
    rDiff(t) = sum(Counts(t,rightpool)) - sum(Counts(t,leftpool));
end
for c = 1:length(poscohs)
    I = abs(coh)==poscohs(c);
    rDiffMean(c) = mean(abs(rDiff(I)));
    rDiffVar(c) = var(rDiff(I));
end
figure(6); subplot(1,2,1); set(gcf,'position',[500 500 1000 420]);
plot(poscohs,rDiffMean,'o-'); title('mean of difference variable (momentary evidence)'); xlabel('unsigned coh');
subplot(1,2,2); plot(poscohs,rDiffVar,'v-'); title('variance of difference variable (momentary evidence)'); xlabel('unsigned coh');
end


%% Back to Mazurek for spike train generation (inverse probability method):

% construct spike density function for each trial according to known time
% course of MT responses, convert to cumulative probability (normalized),
% and use to map uniform random deviates (one for each spike in Counts)
% onto spike times for every trial

% an interesting exercise (which I've yet to do) is check whether this
% method generates the expected cross-correlograms (see e.g. Bair et
% al. J Neurosci 2001). Could be that the 'dummy spike train' from 
% Mazurek is needed for this level of realism.


% baseline = 10; % can be specified independently of tuning curves, but should
               % be higher than response to high-coh motion in antipreferred
               % (null) direction;
       % for now this is set above based on tuning{}               
latency = 50; % ms
risetime = 20; % ms
dT = 1;
tAxis = 1:dT:max(dur)+3*latency;
               % it's 3x because 1x at the beginning, 2x at the end (ramp-down is slower than rise time)

R = nan(nTrials, nNeurons, length(tAxis));
disp('Generating spike trains, may take tens of seconds to minutes...');
tic



for t = 1:nTrials
    for n = 1:nNeurons
                
        meanrate = Rates(t,n);

        sdf=[]; % spike density function
        
            % lay down some baseline rate, from t=1 to vis latency
        sdf(1:latency) = baseline;
            % quick linear rise to mean
        sdf(end+1 : latency+risetime) = linspace(baseline,meanrate,risetime);
%             % steady state with some decay*, like real MT
            decay = 1+0.2*sign(baseline-meanrate); % *up or down depending on whether mean rate is above or below baseline
        sdf(end+1 : latency+dur(t)) = linspace(meanrate, meanrate*decay, dur(t)-risetime);
              % include die-down over 2x latency
        sdf(end+1 : end+2*latency) = linspace(meanrate*decay, baseline, 2*latency);
        
% %         % check what it looks like, if you want
% %         figure; plot(sdf); % CAREFUL! don't run entire loop when uncommented!
                
        % The Mazurek method (I think) assumes zero noise in the driving
        % input (underlying 'rate' or mean of poisson process, see
        % Churchland et al. 2011 Neuron paper "Variance as a Signature...")
        % which is unrealistic even for MT. At the very least it is
        % a function of motion energy, plus any noise added between retina
        % and MT. So adding some ad hoc variability is easily justified
        sdfNoiseSD = 0.2 * sdf; % arbitrary, but scales with spike rate
%         sdfNoiseSD = 0; % arbitrary, but scales with spike rate
        
            % make fluctuations on the time scale of frames (motion energy)?
        % SDFnoise(1:8:dur(t)+2*latency) = SDFnoiseSD*randn(1,length(1:8:dur(t)+2*latency));
        
            % no, for now just every ms; could use actual trials' ME later
        sdf = sdf + sdfNoiseSD.*randn(1,length(sdf));
        sdf(sdf<0)=0;
        
% %         % may also wish to inspect the noisy version
% %         figure; plot(sdf); % CAREFUL! don't run entire loop when uncommented!
                            
        csdf = cumsum(sdf); % cumulative spike density
        cspf = csdf / csdf(end); % cumulative spike probability
            % now we simply draw nspikes uniform random numbers on [0 1] 
            % and use the CSPF as a lookup table for where to place them in
            % time*, redrawing any that fall within a refractory period
            % [*see note below]

            % interpolate to ensure intersect
            % cspfI = interp1(1:length(cspf),cspf,1:0.1:length(cspf));  ??? don't remember what this was about
        
        R(t,n,1:length(sdf)) = 0; % initialize spike train
        % ^ trial, neuron, time    

        sptimes = nan(1,Counts(t,n));
        for s = 1:Counts(t,n) % is there a faster way than spike by spike?
            rnd = rand;
            sp = find(abs(cspf-rnd)==min(abs(cspf-rnd)));
            sp = sp(1); % avoid non-unique times
            if s>1
                while any(abs(sptimes-sp)<2)
                % ^ enforce refractory period by redrawing rnd
                    rnd = rand;
                    sp = find(abs(cspf-rnd)==min(abs(cspf-rnd))); % 
                    sp = sp(1); % avoid non-unique times
                end
            end
            sptimes(s) = sp;
        end
        R(t,n,sptimes) = 1;
        
    end
end
disp('done.');
toc


% check that we achieved the intended spike counts
if plotMT
    Counts2 = nansum(R,3);
    figure(7); plot(Counts(:),Counts2(:),'o',[0 max(max(Counts))+5],[0 max(max(Counts))+5],'k--');
    axis square;
    title('intended vs actual spike counts (should be identical)');
end

% * Note that this method distributes the allocated spikes (determined by
% tuning curves and correlations) throughout the trial, which now
% includes a baseline and ramp-down period. In real experiments, tuning
% curves are usually constructed by counting spikes in a window that
% excludes these periods. This has two consequences: (1) spike rates
% in the appropriate time windows of these constructed spike trains will be
% lower than the prescribed values based on tuning{}, and (2) the average
% spike rate during the baseline period (see PSTHs below) will depend on
% the mean rate for that trial, i.e. will depend on coherence, which of
% course would not happen in a real experiment. I don't think this should
% cause any problems with the simulation, but noting it here just in case.



%% some sanity checks:

if plotMT

% PSTHs, do they look reasonable?
cohsForPSTH = [-0.512 -0.128 -0.032 0.032 0.128 0.512];
clr = cool(length(cohsForPSTH));
figure(8); 
for c = 1:length(cohsForPSTH)
    I = coh==cohsForPSTH(c); % R = [trials, neurons, time]
    plot(nanmean(squeeze(nanmean(R(I,rightpool,1:800),2)))*1000,'Color', clr(c,:)); hold on;
end
title('right pool'); xlabel('time from dots onset (ms)'); ylabel('spikes/s');
legend(cellstr(num2str(cohsForPSTH')),'location','north','orientation','horizontal');
figure(9); 
for c = 1:length(cohsForPSTH)
    I = coh==cohsForPSTH(c);          % R = [trials, neurons, time]
    plot(nanmean(squeeze(nanmean(R(I,leftpool,1:800),2)))*1000,'Color', clr(c,:)); hold on;    
end
title('left pool'); xlabel('time from dots onset (ms)'); ylabel('spikes/s');
legend(cellstr(num2str(cohsForPSTH')),'location','north','orientation','horizontal');


% what about raster plots?
figure(10);
C = coh==0.512;
raster = squeeze(R(C,1,1:400));
[I,J] = find(raster>0);
plot(J,I,'.'); xlabel('time from dots onset (ms)'); ylabel('trial'); ylim([0 size(raster,1)]);
title('spike raster: across trials for a given neuron');
 
figure(11);
C = find(coh==0.512);
raster = squeeze(R(C(1),:,1:400));
[I,J] = find(raster>0);
plot(J,I,'.'); xlabel('time from dots onset (ms)'); ylabel('neuron (pref dir)'); ylim([0 nNeurons]);
title('spike raster: across neurons for a given trial');

end
        
% 
% 
% %% save the results
% disp('saving...');
% 
% clearvars -except R tuning nNeurons nTrials prefDirs latency dir dur coh cohs poscohs tAxis
% 
% tic
% save(sprintf('simMT_nNeu=%d_nTr=%d.mat', nNeurons, nTrials),'-v7.3');
% toc
% 
% disp('done.');
% 



