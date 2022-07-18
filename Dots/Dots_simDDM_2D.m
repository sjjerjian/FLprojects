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
maxdur = 3000;

%% params
% these are very different from 1D model because of the way images_dtb is
% written (sigma is smaller by a factor of 25 drift is greater by about the
% same factor -- why exactly is unclear)

k = 18; % drift rate coeff (conversion from %coh to units of DV)
B = 0.6; % bound height (0.6)
    mu = k*coh; % mean of momentary evidence (drift rate)
% sigma = 0.0376; % SD of momentary evidence -- a free param here but not in images_dtb! used trial and error to find the correct value
sigma = 1; % It's 1 you idiot (adjusted for sqrt(dt), see below).
           % The choice frequencies match up exactly, but the RTs don't, because of the unabsorbed prob!!!
           % (welp, this is only true for certain ranges of k and B, probably for the same reason) 
           
% ME is a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb; can only be changed with an update to the flux file

% convert correlation to covariance matrix
s = [sigma*sqrt(dT/1000) sigma*sqrt(dT/1000)]; % standard deviaton vector
         %^ because variance is sigma^2*dt, so stdev is sigma*sqrt(dt)
% %old:
% s = [sigma sigma]; % standard deviaton vector
         
V = diag(s)*S*diag(s); % equation for covariance matrix
           
theta = 1.2; % threshold for high bet in units of log odds correct (Kiani & Shadlen 09, 14)
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


%% calculate log odds corr maps using Wolpert MoI method
tic
R.t = 0.001:0.001:maxdur/1000; % seconds
R.Bup = B;
R.drift = unique(abs(mu)); % unsigned drift rates
R.lose_flag = 1;
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig
P =  images_dtb_2d(R);
toc

%% bounded evidence accumulation

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
    ME = mvnrnd(Mu,V,maxdur-1);
% %     % check that it has the desired anticorr:
% %     CC = corrcoef(ME); corrcoef_me = CC(1,2)
    dv = [0 0; cumsum(ME)]; % bivariate normrnd
%     dv_all{n} = dv;
    
    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)

    % decision outcome
    cRT1 = find(dv(1:maxdur,1)>=B, 1);
    cRT2 = find(dv(1:maxdur,2)>=B, 1);
    
    % the possibilities are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound,
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = (maxdur)/1000;
        
            % which DV matters for confidence if neither hits bound? 
            % SJ 07/2020 logOddsCorrMap is fixed, so just shift finalV up
            % so that 'winner' did hit bound,
        whichWon = dv(maxdur,:)==max(dv(maxdur,:));
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);
        % ^  shifting the losing dv up by whatever the
        % difference is between the bound and the winning dv

        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
    
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
            
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
            conf(n) = 1 - RT(n) ./ range(P.t);
            pdw(n) = RT(n) < theta;
    end
   
    if isnan(conf(n)), conf(n)=0; end % if dvs are almost overlapping, force conf to zero as it can sometimes come out as NaN

end
toc

% adjust wager probability for base rate of low bets, as seen in data
% ('compresses' the curve, because P(high) varies with coh; see below)
pdw(pdw==1 & rand(length(pdw),1)<alpha) = 0;


%%

% add non-decision time (truncated normal dist)
Tnd = zeros(ntrials,1);
for n = 1:ntrials
    % simple trick for truncating, requires integers (ms)
    while Tnd(n)<=TndMin || Tnd(n)>=TndMax
        Tnd(n) = round(normrnd(TndMean,TndSD));
    end
end
DT = RT; % rename this 'decision time'
RT = DT+Tnd/1000; % convert to s

% quick sanity check that params are reasonable
pCorrect_total = sum(sign(choice)==sign(coh)) / ntrials


%% format data like the real expt

coh(coh==0) = sign(randn)*eps; % should have no actual zeros, but if so, sign them randomly;
                               % this is just to assign a direction and correct/error
data.correct = choice==sign(coh);
data.direction = nan(ntrials,1);
data.direction(coh>0) = 0;
data.direction(coh<0) = 180;
coh(abs(coh)<1e-6) = 0; % now go back to one 'zero'
cohs = unique(coh); % for plotting
data.coherence = abs(coh);
data.scoh = coh;

data.choice = choice;
data.choice(data.choice==-1) = 0; % code elsewhere assumes 0s and 1s
data.RT = RT;
data.PDW = pdw;
data.conf = conf;

conftask = 2; % pdw (2) for now, even though conf rating can be generated here
RTtask = 1; RTCorrOnly = 0;
parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly);

% % TEMP: shouldn't need this; mechanistically the alpha should try to
% % capture what's really going on, which is a coin flip and would result in
% % the 'compression' not the pure offset.
% 
% % adjust PDW to account for a base rate of low-conf bets, 'offset' method
% % cannot just use pdw(pdw==1 & rand(length(pdw),1)<alpha) = 0; because
% % P(high) is not uniform across cohs. Instead need to scale by that prob,
% % and hence this step requires parsedData:
% if alpha>0
%     for c = 1:length(cohs)
%         I = coh==cohs(c);
%         alpha_adj = alpha/parsedData.pHigh(c);
%         pdw(I & pdw==1 & rand(length(pdw),1)<alpha_adj) = 0;
%     end    
%     data.PDW = pdw;
%     parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly);
% end

% plot
wFit = 0; forTalk = 0;
Dots_plot(parsedData,cohs,conftask,RTtask,wFit,forTalk)



%% temp: compare sim to predicted results from images_dtb

% sim:
Pcor = mean([parsedData.pCorrect(6:end) flipud(parsedData.pCorrect(1:6))],2);
RTmean = mean([parsedData.RTmean(6:end)-Tnd(1)/1000 flipud(parsedData.RTmean(1:6))-Tnd(1)/1000],2);

% images_dtb:
figure; set(gcf,'Color',[1 1 1],'Position',[600 600 450 700],'PaperPositionMode','auto'); clf;
suptitle(num2str(sigma));
subplot(2,1,1); plot(R.drift,P.up.p./(P.up.p+P.lo.p),'bo-'); title('Relative prob of corr vs. incorr bound crossed (Pcorr)'); % actually pRight with signed drift
hold on; plot(R.drift,Pcor,'r--x'); % nice. this now matches well with sigma=1
subplot(2,1,2); plot(R.drift,P.up.mean_t,'bo-'); title('mean RT'); xlabel('drift rate');
hold on; plot(R.drift,RTmean,'r--x'); 

% for RT, I assume that mean_t is the mean RT *of trials that crossed the
% bound*, but the empirical (at least sim) data include trials that did not
% hit bound, so adjust model prediction for those:
for c = 1:length(cohs)
    I = data.scoh==cohs(c);
    parsedData.RTmean(c) = mean(RT(I & hitBound')); % redefine to include only hit bound trials
end    
RTmean_new = mean([parsedData.RTmean(6:end)-Tnd(1)/1000 flipud(parsedData.RTmean(1:6))-Tnd(1)/1000],2);
subplot(2,1,2); plot(R.drift,P.up.mean_t,'o-'); title('mean RT'); xlabel('drift rate');
hold on; plot(R.drift,RTmean_new,'g-v'); % close! not quite...


% let's try adjusting the model probs instead:
% probabilty that either bound was crossed:
pHB = P.up.p+P.lo.p;
% under the above hypothesis, the "real" mean RT would be an average of the
% model's reported RT and maxdur, weighted by the probability of hit-bound
% (or not)
realRT_model = pHB.*P.up.mean_t + (1-pHB).*repmat(maxdur/1000,size(P.up.mean_t));

figure; set(gcf,'Color',[1 1 1],'Position',[600 600 450 700],'PaperPositionMode','auto'); clf;
suptitle(num2str(sigma));
subplot(2,1,1); plot(R.drift,P.up.p./(P.up.p+P.lo.p),'bo-'); title('Relative prob of corr vs. incorr bound crossed (Pcorr)'); % actually pRight with signed drift
hold on; plot(R.drift,Pcor,'r--x'); % nice. this now matches well with sigma=1
subplot(2,1,2); plot(R.drift,realRT_model,'bo-'); title('mean RT'); xlabel('drift rate');
hold on; plot(R.drift,RTmean,'r--x'); % also very close! I think we can move on.





