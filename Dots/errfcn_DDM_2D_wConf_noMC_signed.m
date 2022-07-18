function [err,fit,parsedFit] = errfcn_DDM_2D_wConf_noMC_signed(param, guess, fixed, data, options)

% updated model calculations to Steven's method: 
% SJ 10-11-2021 no Monte Carlo simulation for fitting, it's redundant!
% just use model predictions directly

% uses SIGNED cohs for MOI calculation

tic

global call_num

% retrieve the full parameter set given the adjustable and fixed parameters 
param = getParam(param, guess, fixed);

maxdur = 2; % max stimulus duration (s)

k = param(1);
B = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3)); %or negative thetas
alpha = param(4); % base rate of low-conf choices
Tnd = param(5); % fixed Tnd, can't see any other way in this model

% % % sigma = 0.1; % make this a free param? RK didn't need to
% it's not a param at all! hard coded in images_dtb

% use method of images to calculate PDFs of DV and mapping to log odds corr
R.t = 0.001:0.001:maxdur;
R.Bup = B/1.711; % kluge, super weird
R.drift = k * unique(data.scoh); % signed drift rates
R.lose_flag = 1; % we always need the losing densities
R.plotflag = options.plot; % 1 = plot, 2 = plot nicer and export_fig (eg for talk)
R.plotflag = 1; % manual override
P = images_dtb_2d(R); % /WolpertMOI


%% calculate likelihood of the observed data given the model parameters

usetrs_data  = data.correct | data.coherence<1e-6; % use only correct (OR ~0 hdg) trials for RT fits

% trialwise vars
pRight_model_trialwise = nan(length(data.coherence),1);
pHigh_model_trialwise  = nan(length(data.coherence),1);
conf_model_trialwise = nan(length(data.coherence),1);
sigmaConf_data_trialwise = nan(length(data.coherence),1);
RT_model_trialwise = nan(length(data.coherence),1);
sigmaRT_data_trialwise = nan(length(data.coherence),1);

% coh-wise vars (for parsedFit, and some likelihood calculations)
cohs = unique(data.scoh);
n = nan(length(cohs),1);
pRight_model = n;   pHigh_model = n;
meanRT_model = n;   meanRT_data = n;   sigmaRT_data = n;
meanConf_model = n; meanConf_data = n; sigmaConf_data = n;
nCor = n;

for c = 1:length(cohs)
    
    Jdata = data.scoh==cohs(c);

    % pRight_model stores the predictions made by the model for the
    % probability of a rightward choice on each trial condition in the
    % data, and then assigned to the actual data trials of each
    % condition in pRight_model for log likelihood calculation

    % CHOICE
    pRight_model(c) = P.up.p(c)/(P.up.p(c)+P.lo.p(c));        
            % this seems to overestimate the sensitivity (versus sims)
            % could it be because P.up.p is the sum of the bound crossing
            % probability over the full 2 seconds (cdf(end)), whereas some
            % of that probability should 'leak' out as the incorrect bound
            % is crossed? Maybe not, because that prob is given by P.lo.p.
            % E.g. at highest coh, P.lo.p is zero, *also across the 2 sec*,
            % so the P(right) is legitimately 1 in that case. The same
            % logic applies at lower cohs.
    pRight_model_trialwise(Jdata) = pRight_model(c); % copy to trials

    % RT
    nCor(c) = sum(Jdata & usetrs_data);
    if options.RTtask            
        meanRT_model(c) = P.up.mean_t(c) + Tnd; % I think weighting by P.up.p is unnecessary because mean_t already takes it into account            
        RT_model_trialwise(Jdata) = meanRT_model(c); % copy the mean to each trial, for 'fit' struct
        meanRT_data(c) = mean(data.RT(Jdata & usetrs_data));
        sigmaRT_data(c) = std(data.RT(Jdata & usetrs_data)) / sqrt(nCor(c));
        sigmaRT_data_trialwise(Jdata) = sigmaRT_data(c); % copy the SD to each trial, for alternate LL calculation below
    end

    % CONF
    if options.conftask == 1 % sacc endpoint
        warning('code not ready'); keyboard;
        meanConf_model(c) = 2*mean(logistic(confSamples))-1;
        conf_model_trialwise(Jdata & usetrs_data) = meanConf_model(c);
        meanConf_data(c) = mean(data.conf(Jdata & usetrs_data));
        sigmaConf_data(c) = std(data.conf(Jdata & usetrs_data)) / sqrt(nCor(c));
        sigmaConf_data_trialwise(Jdata) = sigmaConf_data(c);
    elseif options.conftask == 2 % PDW        

        Pxt = squeeze(P.up.distr_loser(c,:,:))' .* P.up.p(c);
        Pxt2= squeeze(P.lo.distr_loser(c,:,:))' .* P.lo.p(c);            
        pHigh = sum(Pxt.*(P.logOddsCorrMap>theta)) + sum(Pxt2.*(P.logOddsCorrMap>theta));
        pHigh_model(c) = sum(pHigh);
        pHigh_model_trialwise(Jdata) = pHigh_model(c); % copy to trials
        
        % no, probably need to split by dir in order to get PrightHigh/Low
        
% % %             pRightHigh = cumsum(Ptb(:,2).*(bet_high_tb(:,2)==1)) + ...
% % %                         sum(Pxt(xmesh>0,:).*(bet_high_xt(xmesh>0,:)==1),1)' + ...
% % %                         0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==1))'; % half the probability of dv=0 (assuming a 50/50 guess when that happens)
% % %             pRightLow = cumsum(Ptb(:,2).*(bet_high_tb(:,2)==0)) + ...
% % %                         sum(Pxt(xmesh>0,:).*(bet_high_xt(xmesh>0,:)==0),1)' + ...
% % %                         0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==0))';
% % %             pLeftHigh = cumsum(Ptb(:,1).*(bet_high_tb(:,1)==1)) + ...
% % %                         sum(Pxt(xmesh<0,:).*(bet_high_xt(xmesh<0,:)==1),1)' + ...
% % %                         0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==1))';
% % %             pLeftLow =  cumsum(Ptb(:,1).*(bet_high_tb(:,1)==0)) + ...
% % %                         sum(Pxt(xmesh<0,:).*(bet_high_xt(xmesh<0,:)==0),1)' + ...
% % %                         0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==0))';
        
        % PrightHigh at time t is the the integral of the losing
        % accumulator's density above the theta criterion up to time t,
        % weighted by the cumulative probability that you've crossed the
        % Correct bound by t
        Pxt = squeeze(P.up.distr_loser(c,:,:))'; % density of losing accumulator for this condition
        pRightHigh = P.up.cdf_t(c,:) .* sum(Pxt.*(P.logOddsCorrMap>=theta)); % bound crossing CDF * sum of Pxt above theta
        pRightLow = P.up.cdf_t(c,:) .* sum(Pxt.*(P.logOddsCorrMap<theta));  % bound crossing CDF * sum of Pxt below theta
        Pxt = squeeze(P.lo.distr_loser(c,:,:))';
        pLeftHigh = P.lo.cdf_t(c,:) .* sum(Pxt.*(P.logOddsCorrMap>=theta));
        pLeftLow = P.lo.cdf_t(c,:) .* sum(Pxt.*(P.logOddsCorrMap<theta));
        
        figure;
        subplot(2,2,2); plot(pRightHigh); title('pRightHigh');
        subplot(2,2,4); plot(pRightLow); title('pRightLow');
        subplot(2,2,1); plot(pLeftHigh); title('pLeftHigh');
        subplot(2,2,3); plot(pLeftLow); title('pLeftLow');

        pRightHigh_model(c) = sum(pRightHigh); % marginalize over time [OR use each cond's mean RT??]
        pRightLow_model(c) = sum(pRightLow); 
        pLeftHigh_model(c) = sum(pLeftHigh);
        pLeftLow_model(c) = sum(pLeftLow); 
        
        pRightHigh_model_trialwise(Jdata) = pRightHigh_model(c); % copy to trials
        pRightLow_model_trialwise(Jdata) = pRightLow_model(c);
        pLeftHigh_model_trialwise(Jdata) = pLeftHigh_model(c);
        pLeftLow_model_trialwise(Jdata) = pLeftLow_model(c);

%         if cohs(c)<0
%             Pxt = squeeze(P.lo.distr_loser(c,:,:))'; % density of DV for this condition
%             pRightHigh = sum(Pxt.*(P.logOddsCorrMapL>=theta)); % sum of density above theta [but still is a function of time!]
%             pRightLow = sum(Pxt.*(P.logOddsCorrMapL<theta)); % sum of density above theta [but still is a function of time!]
%         elseif cohs(c)>0
%             Pxt = squeeze(P.up.distr_loser(c,:,:))'; % density of DV for this condition
%             pRightHigh = sum(Pxt.*(P.logOddsCorrMapR>=theta)); % sum of density above theta [but still is a function of time!]
%             pRightLow = sum(Pxt.*(P.logOddsCorrMapR<theta)); % sum of density above theta [but still is a function of time!]
%         else % ?
%         end  
    end
       
end
    
figure; plot(cohs,pRightHigh_model,'b',cohs,pRightLow_model,'r');
% this could be (almost) okay: remember pRightHigh(c) vs pRightLow(c) is
% NOT what we show in the slope effect plots. It's Pright conditioned on
% High vs. Low Bet; not the same thing. For that you need...?



% adjust the probabilities for the base rate of low-conf bets
pHigh_model = pHigh_model - alpha;
pHigh_model(pHigh_model<0) = 0;
pHigh_model_trialwise = pHigh_model_trialwise - alpha;
pHigh_model_trialwise(pHigh_model_trialwise<0) = 0;
% 1D ver does it this way instead (not sure why, but try it if above fails to fit well):
    % pHigh is scaled down in proportion to how close it is 1
    % % pHigh_model = pHigh_model - alpha*pHigh_model;
    % % pHigh_model_trialwise = pHigh_model_trialwise - alpha*pHigh_model_trialwise;

% copy trial params to data struct 'fit', then replace data with model vals
fit = data;
fit = rmfield(fit,'choice');
try fit = rmfield(fit,'correct'); end %#ok<TRYNC>
fit.pRight = pRight_model_trialwise;
if options.conftask==1 % SEP
    fit.conf = conf_model_trialwise;
    fit.PDW = nan(size(fit.conf));
elseif options.conftask==2 % PDW
    fit.conf = pHigh_model_trialwise;
    fit.PDW = pHigh_model_trialwise; % legacy
end
if options.RTtask            
    fit.RT = RT_model_trialwise;
end

% also stored the 'parsed' values, for later plotting
parsedFit = struct();
parsedFit.pRight = pRight_model;
if options.conftask==1 % SEP
    parsedFit.confMean = meanConf_model;
elseif options.conftask==2 % PDW
    parsedFit.pHigh = pHigh_model;
end
if options.RTtask
    parsedFit.RT = meanRT_model;
end



%% Next, until we have bads/ibs_basic working, calculate error using binomial/gaussian assumptions

% convert data vars to logicals
choice = logical(data.choice);
PDW = logical(data.PDW);

% to avoid log(0) issues:
pRight_model_trialwise(pRight_model_trialwise==0) = eps; 
pRight_model_trialwise(pRight_model_trialwise==1) = 1-eps;
pHigh_model_trialwise(pHigh_model_trialwise<=0) = eps; 
pHigh_model_trialwise(pHigh_model_trialwise>=1) = 1-eps;


% CHOICE
% log likelihood of rightward choice on each trial, under binomial assumptions:
LL_choice = sum(log(pRight_model_trialwise(choice))) + sum(log(1-pRight_model_trialwise(~choice)));

% RT
if options.RTtask            
    % likelihood of mean RTs for each coherence, under Gaussian approximation
%     L_RT = 1./(sigmaRT*sqrt(2*pi)) .* exp(-(meanRT_model - meanRT_data).^2 ./ (2*sigmaRT.^2));
    % OR assign means to trials and sum over those, to keep same order of magnitude as other LLs:
    L_RT = 1./(sigmaRT_data_trialwise(usetrs_data)*sqrt(2*pi)) .* exp(-(RT_model_trialwise(usetrs_data) - data.RT(usetrs_data)).^2 ./ (2*sigmaRT_data_trialwise(usetrs_data).^2));

    % ****************
    L_RT(L_RT==0) = min(L_RT(L_RT~=0)); % this seems weird; there shouldn't be any 0 vals anyway...
    % ****************

    LL_RT = nansum(log(L_RT(:))); % sum over all conditions (or trials)
end

% CONF
if options.conftask==1 % SEP
    % likelihood of mean conf ratings for each condition, under Gaussian approximation
    L_conf = 1./(sigmaConf_data*sqrt(2*pi)) .* exp(-(meanConf_model - meanConf_data).^2 ./ (2*sigmaConf_data.^2));
    L_conf(L_conf==0) = min(L_conf(L_conf~=0));
    LL_conf = nansum(log(L_conf(:))); % sum over all conditions
elseif options.conftask==2 % PDW
    % log likelihood of high bet on each trial, under binomial assumptions:
    LL_conf = sum(log(pHigh_model_trialwise(PDW))) + sum(log(1-pHigh_model_trialwise(~PDW)));
end

% total -LL
err = -(LL_choice + LL_conf + LL_RT);
% to add: a method to fit any pair of these and predict the third
% (or any one and predict the other two?)


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
