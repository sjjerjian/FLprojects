function [err,data,fit] = errfcn_DDM_2D_wConf(param, guess, fixed, data, options)

tic

global call_num

% retrieve the full parameter set given the adjustable and fixed parameters 
param = getParam(param, guess, fixed);

ntrials = length(data.choice);
maxdur = 2; % max stimulus duration (s)

k = param(1);
B = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3)); %or negative thetas
alpha = param(4); % base rate of low-conf choices
Tnd = param(5); % need to be a mean + std?
sigma = 0.1; % make this a free param? RK didn't need to

% if Tnd is random var...
% % sdTnd = 60; % fixed SD
% % Tnds = muTnd + randn(ntrials,1).*sdTnd;

% use method of images to calculate PDFs of DV and mapping to log odds corr
R.t = 0.001:0.001:maxdur;
R.Bup = B;
R.drift = k * unique(data.coherence); % takes only unsigned drift rates
R.lose_flag = 1; % we always need the losing densities
R.plotflag = options.plot; % 1 = plot, 2 = plot nicer and export_fig (eg for talk)

P = images_dtb_2d(R); % /WolpertMOI


%% calculate likelihood of the observed data given the model parameters

% *This is ideally done in a similar way as DDM_1D, using the
% probability densities of the terminated DV (marginalized over coherence)
% and log odds corr map. However, because our code originated in the
% context of a continuous confidence judgment and not PDW, we started with
% a simulation-based likelihood calculation (analogous to BADS).

% so, here's a bounded evidence accumulation 'sim' for the 2D model:

choice = nan(ntrials,1);
RT = nan(ntrials,1); 
finalV = nan(ntrials,1);
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
PDW = nan(ntrials,1);

coh = data.scoh; % here we use signed coh

% ME is now a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% this can only be changed with an update to the flux_img file

for n = 1:ntrials
    
    mu = k * coh(n) / 1000; % mean of momentary evidence
    Mu = [mu; -mu]'; % mean vector for 2D DV
    
    % convert correlation to covariance matrix
    V = diag(sigma)*S*diag(sigma);
    
    % DV is cumulative sum of bivariate normal random variates
    dv = [0 0; cumsum(mvnrnd(Mu,V,maxdur*1000))];
    
    % decision outcome
    cRT1 = find(dv(:,1)>=B, 1);
    cRT2 = find(dv(:,2)>=B, 1);
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000 + Tnd; % keep RT in seconds
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound, 
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000 + Tnd;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = 0;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = maxdur + Tnd;
            % this is interesting: which DV matters for confidence
            % if neither hits bound? for now take the (abs) maximum,
            % ie whichever was closest to hitting bound. (alternative
            % would be their average?)
            % SJ 07/2020 finalV is technically distance of loser from bound (when winner
            % hits), so in this case, should also account for where winner is wrt bound 
        whichWon = dv(end,:)==max(dv(end,:));
%         finalV(n) = dv(end,~whichWon); % the not-whichWon is the loser
        % % finalV(n) = mean(dvEnds); 
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);

        hitBound(n) = 0;
        a = [1 0];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000 + Tnd;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 0];
        choice(n) = a(whichWon);
    end
        
    % use map to look up log-odds that the motion is rightward
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
    
    % for post-decision wager, compare LOC to a threshold (theta)
    PDW(n) = logOddsCorr(n) > theta;
end

% output var 'fit' gets the same values as data for the conditions, but the
% simulated trial outcomes for the observables:
fit = data;
fit.choice = choice;
fit.correct = (choice==1 & data.scoh>0) | (choice==0 & data.scoh<0);
fit.RT = RT;
fit.PDW = PDW;


%% Next, until we have bads/ibs_basic working, calculate error using binomial/gaussian assumptions

% convert data vars to logicals
    % (previously the LL calc below simply used choice and PDW from the
    % monte carlo sim! that couldn't be correct.. should be the data :)
choice = logical(data.choice);
PDW = logical(data.PDW);

cohs = unique(data.scoh);
n = nan(length(cohs),1);
pRight_model = n;
pHigh_model = n;
RTmean_fit = n; RTmean_data = n; sigmaRT = n;
nCor = n;

pRight_model_trialwise = nan(ntrials,1);
pHigh_model_trialwise = nan(ntrials,1);

for c = 1:length(cohs)
    J = data.scoh==cohs(c);

    n(c) = sum(J);
    nCor(c) = sum(J & data.correct); % use only correct trials for RT

    pRight_model(c) = sum(J & fit.choice==1) / n(c); % 1 is right
    pRight_model_trialwise(J) = pRight_model(c);

    pHigh_model(c) = sum(J & fit.PDW==1) / n(c); % 1 is high
    pHigh_model_trialwise(J) = pHigh_model(c);
    
    RTmean_fit(c) = mean(fit.RT(J & data.correct));
    RTmean_data(c) = mean(data.RT(J & data.correct));
    sigmaRT(c) = std(data.RT(J & data.correct))/sqrt(nCor(c));
end

% adjust the probabilities for the base rate of low-conf bets
pHigh_model_trialwise = pHigh_model_trialwise - alpha;

% 1D ver does it this way (don't remember why!):
%     % adjust the probabilities for the base rate of low-conf bets:
%     % the idea is that Phigh and Plow each get adjusted down/up in
%     % proportion to how close they are to 1 or 0, respectively
%     Phigh = PrightHigh + PleftHigh;
%     Plow = PrightLow + PleftLow;
%     Phigh2 = Phigh - alpha*Phigh;
%     Plow2 = Plow + alpha*(1-Plow);

% kluge to avoid log(0) issues
% Pr_model(Pr_model==0) = min(Pr_model(Pr_model~=0)); 
% Pr_model(Pr_model==1) = max(Pr_model(Pr_model~=1));
pRight_model_trialwise(pRight_model_trialwise==0) = eps; 
pRight_model_trialwise(pRight_model_trialwise==1) = 1-eps;
pHigh_model_trialwise(pHigh_model_trialwise<=0) = eps; 
pHigh_model_trialwise(pHigh_model_trialwise>=1) = 1-eps;

% log likelihood of rightward choice on each trial, under binomial assumptions:
LL_choice = sum(log(pRight_model_trialwise(choice))) + sum(log(1-pRight_model_trialwise(~choice)));

% log likelihood of high bet on each trial, under binomial assumptions:
LL_conf = sum(log(pHigh_model_trialwise(PDW))) + sum(log(1-pHigh_model_trialwise(~PDW)));

% likelihood of mean RTs for each condition, under Gaussian approximation
L_RT = 1./(sigmaRT*sqrt(2*pi)) .* exp(-(RTmean_fit - RTmean_data).^2 ./ (2*sigmaRT.^2));
L_RT(L_RT==0) = min(L_RT(L_RT~=0));
LL_RT = nansum(log(L_RT(:))); % sum over all conditions

% total -LL
err = -(LL_choice + LL_conf + LL_RT);


%% print progress report!
if options.feedback
    fprintf('\n\n\n****************************************\n');
    fprintf('run %d\n', call_num);
    fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\talpha= %g\n\tTnd= %g\n', k, B, theta, alpha, Tnd);
    fprintf('err: %f\n', err);
end
if options.feedback==2 && strcmp(options.fitMethod,'fms')
    figure(options.fh);
    plot(call_num, err, '.','MarkerSize',14);
    drawnow;
end
call_num = call_num + 1;
toc

end


% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam ( param1 , guess , fixed )
  param2(fixed==0) = param1(fixed==0);  %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess
end
