%% %YYY% 07/08/22 Temporally sequential MT_Build

% Each neuron has a predetermined spRate for a given trial, stored in
% 'Rates', generated above. From each trial spRate we will generate a spike
% vector that details when in time each spike occured on the given trial.
% Spike vector is currently generated using cumsum(sdf) as a look up table
% for the time in trial with a spiking probability that best matches a rand
% generated probability value 'rnd'. %YYY%
    % This generates spike count, rastorlike, vectors by nonsequentially jumping around
    % along time axis, forward and backwards. Thus, unable to track
    % evolving decisioni process through time. Must find way to create
    % spike vector linearly moving forward in time.
%How? Given a vector of spiking probability at every time point in
%trial...can loop through time iteratively (in whatever bin width) and ask
%at each timepoint for each neuron, is a spike generated? Spike is
%generated when the probability of spiking matches the rnd generated value!
% The output will be a matrix of all neurons (rows) storing all spikes at a
% given timepoint in a trial (neurons x trial x time)

% The value in generating spikes sequentially across time is to later slip
% in CP readout and "feedback effects" that act to alter the spike
% responses later in trial. Build in pool classification as well... pull
% from MT_test ?


%% build a simulated MT population and generate its spiking activity to a 
% standard set of random-dot motion stimuli
% CF 2014-2020

clear
close all
    
%% build an experiment
    % Define number of trials to simulate and available trial stim coherences (cohs), Each trial is assigned a duration (dur) and stim coherence (coh) value

tic

% plot flags
plotMT = 1;
plotCorr = 0;
plotLOmap = 0; % for DDM plot log-odds corr map or not (set to 2 to see P(DV) also)


% noiseModel = 'indepPoisson';
% noiseModel = 'corrMat_Mazurek02';
noiseModel = 'corrMat_Shadlen96';

nTrials = 250;

% signed motion coherence; negative is leftward
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -1e-8 1e-8 0.032 0.064 0.128 0.256 0.512];  % All possible coherence levels, signed
coh = randsample(cohs,nTrials,'true')'; % Assign a coh value to each of 1000 trials

dir = nan(size(coh));   % Create 'dir' vector 1 column, 1 row for each trial (nTrials = 1000)
dir(coh<0) = 180;   % Index to trials (row) with coh < 0. Assign these trials value 180 = Left is correct answer 
dir(coh>0) = 0;     % Same as above, but now Right is correct answer

% viewing duration drawn from truncated exponential distribution
% dur = zeros(nTrials,1); % Vector to contain the duration value of each trial %% USE if not using truncexpdur function!
% meandur = 360; mindur = 80; maxdur = 900;   % Duration parameters
% dur = truncexpdur(meandur,mindur,maxdur,nTrials); % Creates ms duration for each trial, trial duration values form an exponential distribution
%  truncexpdur does below
% for n = 1:nTrials
%     while dur(n)<=mindur || dur(n)>=maxdur % simple trick for truncating
%         dur(n) = round(exprnd(meandur));
%     end
% end

% OR set dur to 1s (or any fixed value) -- btw makes it easier to verify fano factor (below)
dur(1:nTrials,1) = 650;

toc

        disp('Task structure defined');
%% build MT population 
    % Each neuron assigned a preferred motion direction between 0,360 degrees
    % Figure2: Plot a neurons tuning curve response at all coh levels and all directions of motion
    % Figure3: Plots rightpool and leftpool frRates to al lstim directions at single coh level
    

tic

% nNeurons = 36;
nNeurons = 180;
% nNeurons = 360;

dirAxis = 0:2:359;
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
%
%
% SDF stores the vector of each individual neurons sdf row by row
% For each time (ms) spike is generated or not, simultaneously for all neurons.
  
latency = 50; % ms
risetime = 20; % ms
dT = 1;
tAxis = 1:dT:max(dur)+3*latency; % it's 3x because 1x at the beginning, 2x at the end (ramp-down is slower than rise time)

SDF = nan(nTrials, nNeurons, length(tAxis));  % spProb for every neuron, on each trial, at every time point within that trial

Ry = nan(nTrials, nNeurons, length(tAxis));  % [trial x neuron x time(ms)] %Fill Ry with rastor vectors. Spike or no spike at each time point, of each trial, for each neuron. Leave as nans values past dur(tr) to avoid averaging into pop 
% spTimes = nan(nTrials,nNeurons);   % Counts(trials x neuorns) indexes to the number of spikes counted on trial 'tr', for neuron 'n'.

disp('Generating spike trains, may take tens of seconds to minutes...');

tic  % Took ~50 seconds to complete 180 neurons, 350 trials, 650 trial dur,1:1:ms
for tr = 1:nTrials
    for n = 1:nNeurons

        meanrate = Rates(tr,n);   % Not sure why called "meanrate"... is just spRate for neuron 'n' on trial 'tr'

        sdf=[]; % spike density function

        % lay down some baseline rate, from tr=1 to vis latency
        sdf(1:latency) = baseline;
        % quick linear rise to mean
        sdf(end+1 : latency+risetime) = linspace(baseline,meanrate,risetime);
        %             % steady state with some decay*, like real MT
        decay = 1+0.2*sign(baseline-meanrate); % *up or down depending on whether mean rate is above or below baseline
        sdf(end+1 : latency+dur(tr)) = linspace(meanrate, meanrate*decay, dur(tr)-risetime);
        % include die-down over 2x latency
        sdf(end+1 : end+2*latency) = linspace(meanrate*decay, baseline, 2*latency);

        % %         % check what it looks like, if you want ('sdf' is unique for each trial, each is stored in 'SDF')
        % %         figure; plot(sdf); % CAREFUL! don'tr run entire loop when uncommented!
        %
        %         % Add some NOISE to synthetic neural responses
        sdfNoiseSD = 0.2 * sdf; % arbitrary, but scales with spike rate
        % %             % make fluctuations on the time scale of frames (motion energy)?
        %         SDFnoise(1:8:dur(tr)+2*latency) = SDFnoiseSD*randn(1,length(1:8:dur(tr)+2*latency));
        % %             % no, for now just every ms; could use actual trials' ME later (ME = motion energy?)
        sdfNoise = sdf + sdfNoiseSD.*randn(1,length(sdf));
        sdfNoise(sdfNoise<0)=0;
        %
        % % %         % may also wish to inspect the noisy version
        % % %         figure; plot(sdf); % CAREFUL! don'tr run entire loop when uncommented!
        %
        %         csdf = cumsum(sdf); % cumulative spike density
        %         csdfNoise = cumsum(sdfNoise); % cumulative spike density
        %         cspf = csdf / csdf(end); % cumulative spike probability
        %         cspfNoise = csdfNoise / csdfNoise(end); % cumulative spike probability


        SDF(tr,n,1:length(sdf)) = sdfNoise; % Stores the spProb for each neuron, at each timepoint, on every trial.
        %         Can change value ^ between 'sdf' or 'sdfNoise' as desired
    end

%%%%% SDF built! out of 'n' loop, still in 'tr' loop
    % At each ms in trial tr, ask whether each neuron spiked or not.
    % Probabilistically determined by sdf value for that particular neuron, at time ms.
    % A decision variable 'dv' is generated from the
    % population cumulative spike counts measured from time 1:ms
    % Whether this dv supports R or L choice determines modulatory
    % effects on both R and Lpool probabilistic firing rates ('lambda'), for time ms+1.

    trTime = dur(tr)+3*latency; % Time length of each trial

    for ms = 1:dT:trTime % At each ms timepoint compare 'rnd' to each neurons spike probability. If the values are close in value spike occurs. Define **close**: In simMT_build difference between rnd and cspf value is minimized. Capture this minimized value for each timepoint.
        % if wanna do bigger time steps, must account for cumulative
        % probaility within the time bin size so cumsum(SDF(tr,:,ms:ms+dT-1)
        % Each 1sec has a spProb stored in SDF. CSPF stores the cumulative and
        % normalized (i.e. relative). cumsum() = value 1, value 1 +v2, v1 +v2 +v3,
        % end(value) = 1 when normalized (i.e. divide all by max value, (also the end value)
        % Using the spProb for a single sec we can run that through poissrnd() 1000
        % times (for number of ms in a second), arriving us at a spike vector for that 1 sec
        % OR %
        % We can divide the 1s spProb value by 1000 to scale to ms. Now have a
        % constant spProb for every 1ms within the 1sec timebin. Can do for every
        % 20ms if divide 1s spProb value by 50 (1000/20) rather than 1000. dT = #ms
        % bins 1s is scaled into. ---> This latter description is what we do. See
        % how lambda is defined as SDF / 1000/dT

        % Can try generating lambda from diff vals: Use SDF, CSDF, CSPF. Can also
        % try cumsum(lambda) but needa make sure cumulative summing within same
        % neuron, i.e.
        %
        % tried
        % doing cumsum(lambda) obvi didn't work cuz cumsum function adds together
        % values in same vector, thus combining lambda values across neurons, not
        % across time as desired. Thus, instead updated lambdas value each ms
        % adding to previous times lambda.

        lambda = SDF(tr,:,ms) / (1000/dT);  % Get "OG" lambda value for time 'ms'
                                            % This lambda value will be modified (or not) in if statement below
        % lambda = spProb in specified timebin, of size dT
        % SDF holds spike probability of each neuron for specific trial at
        % specific timepoint. That timepoint is = 1second. but we wanna
        % know all the spikes that occur within that 1 second (dT).
        % Distribute SDF spProb/ 1sec --> spProb 1sec/ dT, if want to know
        % spProb for every ms set dT=1 --> 1sec=1000ms so, SDF / (1000ms/1)

%          % READ HERE for Updating value of lambda   --->> This code is in the if statements below but for legibility I leave it commented here all together (expand the comment)               
%         if ms == 1
%             Ry(tr,:,ms) = poissrnd(lambda);
%         else % lambda is defined by sdf for each neuron at every 'ms'. This pre-defined lambda is modified at the current (ms >1) as a function of the state of DV at (ms-1)
%             lambdaL = lambda +1; %*some operation*, currently "+1"% ; % Old + new what operation modulates lambda appropriately relative to logOdds value
%             lambdaR = lambda +1; %*some operation*,
% 
%             Ry(tr,rightpool,ms) = poissrnd(lambdaR);
%             Ry(tr,leftpool,ms) = poissrnd(lambdaL);  
%         end       

        % Initialize Decode and PDW task variables
        if ms == 1 && tr == 1 %  Only want to run these operations one time total
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

            rightpool_idx = prefDirs(rightpool==1);
            leftpool_idx = prefDirs(leftpool==1);

        elseif ms == 1 % Run once for each trial
            Ry(tr,:,ms) = poissrnd(lambda); % Generate spike or no spike for time ms ==1

            rightspikes(dur(tr)+latency+1:end) = NaN; % Start vector after latency then fill times beyond RT with Nans
            leftspikes(dur(tr)+latency+1:end) = NaN;
        else % i.e., for any trial and all ms>1hh
            if RgreaterL ==1
                lambdaLoser = lambda; % placeholder operation; % Old + new what operation modulates lambda appropriately relative to logOdds value
                lambdaWinner = lambda; % placeholder operation

                lambda(rightpool) = lambdaWinner(rightpool); % lambda is a vector of 360 values
                lambda(leftpool) = lambdaLoser(leftpool);
                Ry(tr,:,ms) = poissonrnd(lambda);
            else % i.e., R is not greater L, i.e., LgreaterR ==1
                lambdaLoser = lambda; % *probR placeholder operation; % Old + new what operation modulates lambda appropriately relative to logOdds value
                lambdaWinner = lambda; % *probL placeholder operation

                lambda(leftpool) = lambdaWinner(leftpool); 
                lambda(rightpool) = lambdaLoser(rightpool);
                Ry(tr,:,ms) = poissrnd(lambda);
            end
        end

        rightspikes(tr,ms) = sum(sum(Ry(tr,rightpool,1:ms),3,'omitnan')); % Current spike count total for each pool, counts through till time ms
        leftspikes(tr,ms) = sum(sum(Ry(tr,leftpool,1:ms),3,'omitnan')); % Try without 'omitnan' so length matches up

        poolingNoise = round(sigmaPooling*randn(1,size(rightspikes,2)));

        diff = rightspikes(tr,:) - leftspikes(tr,:) + poolingNoise;
        dv = cumsum(diff); % cumulative difference between pool spCounts
        % want to look up in logOddsMapL/R the value that corresponds to a
        % particular dv (the dv closest to our current dv).
        % Find index to the closest vAxis value to current dv.
        % Enter this index and current ms as coordinates to locate prob value from
        % logOddsMapR/L.

        dv_lookup = find(min(abs(vAxis-dv(ms)))); % Find index to vAxis value closest to current dv (at time 'ms').
        probL = abs(logOddsMapL(dv_lookup,ms)); % Index to corresponding likelihood value for L and R
        probR = abs(logOddsMapR(dv_lookup,ms)); % Compare these values to see which is more likely
        RgreaterL = probR > probL;  % Yes or No?
    end

    % Sim PDW task behavior, now that a full trial spike vector is available

        % RT
        cRT = find(abs(dv)>=B, 1, 'first');
        if isempty(cRT)
            RT(tr) = dur(tr);
            finalV(tr) = dv(RT(tr));
            hitBound(tr) = 0;
        else
            RT(tr) = cRT;
            finalV(tr) = B*sign(dv(RT(tr))); % cap bound crossings at bound value
            hitBound(tr) = 1;
        end

        % choice
        if finalV(tr)==0
            choice(tr) = sign(rand-0.5);
        else
            choice(tr) = sign(finalV(tr));
        end

        % look up the expected log odds corr for `14567this trial's V and T
        whichV = find(abs(vAxis-finalV(tr))==min(abs(vAxis-finalV(tr)))); % Finds vAxis value closest to finalV. %YYY% seems redundant..? Why not just find(min(abs(vAxis-V(tr))))
        whichT = find(abs(tAxis-RT(tr))==min(abs(tAxis-RT(tr))));
        logOddsCorr(tr) = logOddsCorrMap(whichV(1), whichT(1));

        % post-decision wager:
        if logOddsCorr(tr) >= theta % bet high
            pdw(tr) = 1;
        elseif logOddsCorr(tr) < theta % bet low
            pdw(tr) = 0;
        else
            keyboard
        end

        choice(tr) = (choice(tr)+1)/2; % convert to 0/1

        %    parfor_progress;
end
toc
disp('Spike vectors generated');

% for PDW RT behaviour, Not sure if necessary? Maybe for behavioral plots..
TndMean = 300;
TndSD = 50;
TndMin = TndMean/2;
TndMax = TndMean+TndMin;
Tnd = zeros(nTrials,1);
for tr = 1:nTrials
    while Tnd(tr)<=TndMin || Tnd(tr)>=TndMax % simple trick for truncating
        Tnd(tr) = round(normrnd(TndMean,TndSD));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd;

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



%% Assesment plots of Neural spike trains

tic
if plotMT
%     SDFCounts = sum(Ry,3,'omitnan');  % Good counts but PSTH plots falling activity
%   LamCounts = sum(Ry,3,'omitnan');  % Interesting.. Compensates well for falling PSTH from including 0's in average rather than Nans.

    Counts2 = sum(Ry,3,'omitnan');
    figure(7); plot(Counts(1:tr-1,:),Counts2(1:tr-1,:),'o',[0 max(max(Counts))+5],[0 max(max(Counts))+5],'k--');
    axis square;
    title('intended vs actual spike counts (should be identical)');

    % PSTHs, do they look reasonable?
    cohsForPSTH = [-0.512 -0.128 -0.032 0.032 0.128 0.512];
    clr = cool(length(cohsForPSTH));
    figure(8); 
    for c = 1:length(cohsForPSTH)
        I = coh(1:tr-1) == cohsForPSTH(c); % Ry = [trials, neurons, time]
        plot(mean(squeeze(mean(Ry(I,rightpool,1:800),2,'omitnan')),'omitnan')*1000,'Color', clr(c,:)); hold on;
    end
    title('right pool'); xlabel('time from dots onset (ms)'); ylabel('spikes/s');
    legend(cellstr(num2str(cohsForPSTH')),'location','north','orientation','horizontal');
    figure(9); 
    for c = 1:length(cohsForPSTH)
        I = coh(1:tr-1) == cohsForPSTH(c); % Ry = [trials, neurons, time]
        plot(mean(squeeze(mean(Ry(I,leftpool,1:800),2,'omitnan')),'omitnan')*1000,'Color', clr(c,:)); hold on;    
    end
    title('left pool'); xlabel('time from dots onset (ms)'); ylabel('spikes/s');
    legend(cellstr(num2str(cohsForPSTH')),'location','north','orientation','horizontal');

    % raster plots
    figure(10);
    C = coh==0.512;
    raster = squeeze(Ry(C,1,1:400));
    [I,J] = find(raster>0);
    plot(J,I,'.'); xlabel('time from dots onset (ms)'); ylabel('trial'); ylim([0 size(raster,1)]);
    title('spike raster: across trials for a given neuron');
     
    figure(11);
    C = find(coh==0.512);
    raster = squeeze(Ry(C(1),:,1:400));
    [I,J] = find(raster>0);
    plot(J,I,'.'); xlabel('time from dots onset (ms)'); ylabel('neuron (pref dir)'); ylim([0 nNeurons]);
    title('spike raster: across neurons for a given trial');
end
toc
disp('Neural assesments plotted');

%% save the results
% For use in simMT_test
% Make sure to initiate MT_test with nNeuron and nTrial values matching MT_build

% disp('saving...');
% % 
% clearvars -except tr ms lambda lambdaWinner lambdaLoser finalV hitBound ogOddsMapL logOddsMapR logOddsCorrMap Mu Sigma pdw choice rightpool leftpool Ry tuning nNeurons nTrials prefDirs latency dir dur coh cohs poscohs tAxis
% % 
% tic
% save(sprintf('simMT_nNeu=%d_nTr=%d.mat', nNeurons, nTrials),'-v7.3');
% toc
% % 
% disp('file saved.');

%% 


%% sanity check: choice probability (Britten et al. 1996)
% clear;
% load('simMT_nNeu=360_nTr=1000.mat')
I = abs(coh)<0.0001;
chooseright = I & choice==1; sum(chooseright);
chooseleft = I & choice==0; sum(chooseleft);
% first try based on full spike count
spCount = squeeze(sum(Ry,3,'omitnan'));
k=1;
for g = find(rightpool)
    CP(k,1) = rocN(spCount(chooseright,g), spCount(chooseleft,g), 50);
    k = k+1;
end
meanCP = mean(CP)

%% Choice Probability/ Decision readouts from Spike vectors!
% clear;
% load('simMT_nNeu=360_nTr=1000.mat')

CP_high = nan(size(1:nNeurons));

CP_low = nan(size(1:nNeurons));

CP_all = nan(size(1:nNeurons));

maxCoh = .0001;
right = 1; left = 0;
Counts2 = sum(Ry,3,'omitnan');

tic
for tr = 1: nTrials
    for n = 1:nNeurons
        if rightpool(n) == 1
            pref = right;
        else
            pref = left;
        end
        
        Rates2 = Counts2(tr,:) ./(dur/1000); % Stores 1 spRate value for each neuron on every trial, i.e. in each column of each row.
        cohIdx = abs(coh) < maxCoh;
       
        prefHi = cohIdx & choice ==pref & pdw == 1; % PDW Hi; Logical array indexing to trials of prefHi choice outcome
            X = Rates2(prefHi,n); % Index to SpRates of neuron 'n' on trials where choice = pref and pdw = Hi
        nullHi = cohIdx & choice ~=pref & pdw == 1; % PDW Hi;
            Y = Rates2(nullHi,n); % null-high
        CP_high(n) = rocN(X,Y,100);
    
    
        prefLo = cohIdx & choice ==pref & pdw == 0; % PDW Lo;
          X = Rates2(prefLo,n); % pref-low
        nullLo = cohIdx & choice ~=pref & pdw == 0; % PDW Lo;
          Y = Rates2(nullLo,n); % null-low
        CP_low(n) = rocN(X,Y,100);
    
        preff = cohIdx & choice ==pref; % No PDW;
          X = Rates2(preff,n); % pref-all
        null = cohIdx & choice ~=pref; % No PDW;
          Y = Rates2(null,n); % null-all
        CP_all(n) = rocN(X,Y,100);
    
    end
end
toc

% %     % Can also do ConfP
