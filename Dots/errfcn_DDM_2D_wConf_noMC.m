function [err,fit,parsedFit] = errfcn_DDM_2D_wConf_noMC(param, guess, fixed, data, options)

% updated model calculations to Steven's method: 
% SJ 10-11-2021 no Monte Carlo simulation for fitting, it's redundant!
% just use model predictions directly

% uses UNSIGNED cohs for MOI calculation

tic

global call_num

% retrieve the full parameter set given the adjustable and fixed parameters 
param = getParam(param, guess, fixed);
max_dur = max(data.dur)/1000; % max stimulus duration (s)

k = param(1);
B = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3)); % or negative thetas
alpha = param(4); % base rate of low-conf choices
Tnd = param(5); % fixed Tnd, can't see any other way in this model

% use method of images to calculate PDFs of DV and mapping to log odds corr
R.t = 0.001:0.001:max_dur;
R.Bup = B;
R.drift = k * unique(data.coherence); % takes only unsigned drift rates
R.drift_freq = hist(data.coherence,unique(data.coherence))'/length(data.coherence); % relative frequency of each coh(drift)
R.lose_flag = 1; % we always need the losing densities
R.plotflag = options.plot; % 1 = plot, 2 = plot nicer and export_fig (eg for talk)
%     R.plotflag = 0; % manual override
P = images_dtb_2d(R); % /WolpertMOI


%% calculate likelihood of the observed data given the model parameters

% usetrs_data = data.correct | data.coherence<1e-6; % use only correct (OR ~0 hdg) trials for RT fits
usetrs_data = true(size(data.correct)); % try all trials

% coh-wise vars (for parsedFit, and some likelihood calculations)
cohs = unique(data.scoh);
n = nan(length(cohs),1);
pRight_model = n;
pHigh_model = n;
pRightHigh_model = n;
pRightLow_model = n;
pHighCorr_model = n;
pHighErr_model = n;
meanRT_model = n;
meanRThigh_model = n;
meanRTlow_model = n;
meanRT_data = n;
meanRThigh_data = n;
meanRTlow_data = n;
sigmaRT_data = n;
sigmaRThigh_data = n;
sigmaRTlow_data = n;
meanConf_model = n;
meanConf_data = n;
sigmaConf_data = n;
nCor = n;

% trialwise vars
m = nan(length(data.coherence),1);
pRight_model_trialwise = m;
pHigh_model_trialwise  = m;
pRightHigh_model_trialwise = m;
pRightLow_model_trialwise = m;
pLeftHigh_model_trialwise = m;
pLeftLow_model_trialwise = m;
conf_model_trialwise = m;
RT_model_trialwise = m;
RThigh_model_trialwise = m;
RTlow_model_trialwise = m;
sigmaRT_data_trialwise = m;
sigmaRThigh_data_trialwise = m;
sigmaRTlow_data_trialwise = m;

for c = 1:length(cohs) % loop through signed cohs, because that is what defines a trial in the data
    
    % index for images_dtb (unsigned coh/drift corresponding to this signed coh)                    
    uc = abs(cohs(c))==(cohs(cohs>=0));  % *** marking differences in signed vs. unsigned ver ***
      
    % grab the distributions from images_dtb:
    Pxt1 = squeeze(P.up.distr_loser(uc,:,:))'; % losing race pdf given correct
    Pxt2 = squeeze(P.lo.distr_loser(uc,:,:))'; % losing race pdf given incorrect

    % CHOICE & PDW
    if cohs(c)<0 % leftward motion
        % if coh is negative, then pRight is p(incorrect), aka P.lo
        pRight = P.lo.p(uc)/(P.up.p(uc)+P.lo.p(uc)); % normalized by total bound-crossing probability 
        pLeft = 1-pRight; % (evidently, the unabsorbed probability can be ignored when it comes to pRight)

        % calculate probabilities of the four outcomes: right/left x high/low
        % these are the intersections, e.g. P(R n H)
        pRightHigh = sum(sum(Pxt2.*(P.logOddsCorrMap>=theta)));
        pRightLow = sum(sum(Pxt2.*(P.logOddsCorrMap<theta)));
        pLeftHigh = sum(sum(Pxt1.*(P.logOddsCorrMap>=theta)));
        pLeftLow =  sum(sum(Pxt1.*(P.logOddsCorrMap<theta)));                                                             
    else % rightward motion
        % if coh is positive, then pRight is p(correct), aka P.up
        pRight = P.up.p(uc)/(P.up.p(uc)+P.lo.p(uc)); % normalized by total bound-crossing probability
        pLeft = 1-pRight;

        % calculate probabilities of the four outcomes: right/left x high/low
        % these are the intersections, e.g. P(R n H)
        pRightHigh = sum(sum(Pxt1.*(P.logOddsCorrMap>=theta)));
        pRightLow = sum(sum(Pxt1.*(P.logOddsCorrMap<theta)));
        pLeftHigh = sum(sum(Pxt2.*(P.logOddsCorrMap>=theta)));
        pLeftLow =  sum(sum(Pxt2.*(P.logOddsCorrMap<theta)));                                                     
    end
    
    %ensure total prob sums to one (when it doesn't, it's because of unabsorbed prob)
    % THIS STEP WASN'T THOUGHT TO BE NECESSARY FOR PARAM RECOV, SO BE
    % CAREFUL AND REMOVE IT IF IT MUCKS IT UP
    Ptot = pRightHigh + pRightLow + pLeftHigh + pLeftLow;
    pRightHigh = pRightHigh/Ptot;
    pRightLow = pRightLow/Ptot;
    pLeftHigh = pLeftHigh/Ptot;
    pLeftLow = pLeftLow/Ptot;

    % copy intersections to trials, for multinomial likelihood
    Jdata = data.scoh==cohs(c); % select trials of this coh
    pRightHigh_model_trialwise(Jdata) = pRightHigh;
    pRightLow_model_trialwise(Jdata) = pRightLow;
    pLeftHigh_model_trialwise(Jdata) = pLeftHigh;
    pLeftLow_model_trialwise(Jdata) = pLeftLow;
    
    % by definition of conditional probability: P(A|B) = P(A n B) / P(B)
    pHigh_Right = pRightHigh / pRight;
    pLow_Right = pRightLow / pRight; 
    pHigh_Left = pLeftHigh / pLeft;
    pLow_Left = pLeftLow / pLeft;  
    
    % by law of total probability
    pHigh = pHigh_Right*pRight + pHigh_Left*pLeft;
    pLow = pLow_Right*pRight + pLow_Left*pLeft; % why isn't this 1-pHigh? probability is so weird.

    % by Bayes' rule:
    pRight_High = pHigh_Right * pRight / pHigh;
    pRight_Low = pLow_Right * pRight / pLow;
    pLeft_High = pHigh_Left * pLeft / pHigh;
    pLeft_Low = pLow_Left * pLeft / pLow; %#ok<NASGU>

    % adjust the probabilities for the base rate of low-conf bets (alpha)
    pHigh_wAlpha = pHigh - alpha*pHigh;

    % Lastly, the PDW split; use Bayes to adjust P(H|R) and P(H|L)
    pHigh_Right_wAlpha = pRight_High * pHigh_wAlpha / pRight;
    pHigh_Left_wAlpha = pLeft_High * pHigh_wAlpha / pLeft;
    if cohs(c)<0 % leftward motion
        pHigh_Corr = pHigh_Left_wAlpha;
        pHigh_Err = pHigh_Right_wAlpha;
    else % rightward motion
        pHigh_Corr = pHigh_Right_wAlpha;
        pHigh_Err = pHigh_Left_wAlpha;    
    end
    
    % copy to vectors for parsedFit struct
    pRight_model(c) = pRight;
    pHigh_model(c) = pHigh_wAlpha; % includes the effect of alpha
    pRightHigh_model(c) = pRight_High; % does NOT include the effect of alpha
    pRightLow_model(c) = pRight_Low; % does NOT include the effect of alpha
    
% % %     pLeftHigh_model(c) = pLeftHigh; % MAYBE TEMP: SEE LIKELIHOOD
% % %     pLeftLow_model(c) = pLeftLow;  % MAYBE TEMP: SEE LIKELIHOOD
    
    pHighCorr_model(c) = pHigh_Corr; % includes the effect of alpha
    pHighErr_model(c) = pHigh_Err; % includes the effect of alpha

    % copy to trials for fit struct, and likelihood
    pRight_model_trialwise(Jdata) = pRight_model(c);
    pHigh_model_trialwise(Jdata) = pHigh_model(c);
    pRightHigh_model_trialwise(Jdata) = pRightHigh_model(c);
    pRightLow_model_trialwise(Jdata) = pRightLow_model(c);
% % %     pHighCorr_model_trialwise(Jdata) = pHighCorr_model(c);
% % %     pHighErr_model_trialwise(Jdata) = pHighErr_model(c);
        
    % RT
    nCor(c) = sum(Jdata & (data.correct | data.coherence<1e-6));
    if options.RTtask            
        
        % for param-recov, model RT needs adjusted to account for unabsorbed
        % probability (trials that don't hit the bound, usually only at low
        % cohs, when bound is high, etc). Turns out the calculation of
        % P.up.mean_t only applies to absorbed prob (bound-crossings), so
        % we simply adjust it by averaging mean_t with maxdur, weighted by
        % one minus the probability of bound crossing
            % (this was a cool exercise when it worked, but seems easier to
            % exclude those trials and pretend there is no unabsorbed prob;
            % also more consistent with real behavior where there's always a
            % choice before max_dur is reached  --thanks to MVL for this)
        if options.ignoreUnabs
            pHB = 1;
        else
            pHB = P.up.p(uc)+P.lo.p(uc); % prob of crossing either bound
        end
        
        meanRT_model(c) = pHB*P.up.mean_t(uc) + (1-pHB)*max_dur + Tnd; % weighted average
 
        % for RT conditioned on wager, it seems we don't need to do any
        % scaling by Ptb; the correct distribution is captured by the
        % conditional Pxts, we just sum them, dot product w t and normalize
        PxtAboveTheta = sum(Pxt1.*(P.logOddsCorrMap>=theta)); % shouldn't matter if Pxt1 or 2, it's symmetric
        PxtBelowTheta = sum(Pxt1.*(P.logOddsCorrMap<theta));
        meanRThigh = PxtAboveTheta * R.t' / sum(PxtAboveTheta);
        meanRTlow = PxtBelowTheta * R.t' / sum(PxtBelowTheta);
        
        % BUT these also need to factor in unabsorbed probability
        % (ie weighted average with max_dur). Intuitively, we know that we    
        % must scale pHB up for high bets, down for low bets. But I was
        % unable to calculate the adjusted pHBs and gave up, for now (see
        % old notes)
        % Instead, just set both to pHB, which is close enough when pHB>.98
        pHBhigh = pHB;
        pHBlow = pHB;

        meanRThigh_model(c) = pHBhigh*meanRThigh + (1-pHBhigh)*max_dur + Tnd;
        meanRTlow_model(c) = pHBlow*meanRTlow + (1-pHBlow)*max_dur + Tnd;
        
        % copy to trials
        RT_model_trialwise(Jdata) = meanRT_model(c);
        RThigh_model_trialwise(Jdata) = meanRThigh_model(c);
        RTlow_model_trialwise(Jdata) = meanRTlow_model(c);
        
        % grab the mean and sd(se?) for corresponding data
        I = Jdata & usetrs_data;
        meanRT_data(c) = mean(data.RT(I));
        sigmaRT_data(c) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
        sigmaRT_data_trialwise(Jdata) = sigmaRT_data(c); % copy to trials, for alternate LL calculation below
                                                         % (should this be Jdata or I?)            
        % repeat for high/low
        I = Jdata & usetrs_data & data.PDW_preAlpha==1;
        meanRThigh_data(c) = mean(data.RT(I));
        sigmaRThigh_data(c) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
        sigmaRThigh_data_trialwise(Jdata) = sigmaRThigh_data(c); % copy to trials, for alternate LL calculation below
        I = Jdata & usetrs_data & data.PDW_preAlpha==0;
        meanRTlow_data(c) = mean(data.RT(I));
        sigmaRTlow_data(c) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM?
        sigmaRTlow_data_trialwise(Jdata) = sigmaRTlow_data(c); % copy to trials, for alternate LL calculation below

    end
    
end

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
    parsedFit.pRightHigh = pRightHigh_model;
    parsedFit.pRightLow = pRightLow_model;
    parsedFit.pHighCorr = pHighCorr_model;
    parsedFit.pHighErr = pHighErr_model;
end
if options.RTtask
    parsedFit.RTmean = meanRT_model;
    if options.conftask==2 % PDW
        parsedFit.RTmeanHigh = meanRThigh_model;
        parsedFit.RTmeanLow = meanRTlow_model;
    end
end


%% Next, calculate error (negative log-likelihood) using binomial/gaussian,
% assumptions

% convert data vars to logicals
choice = logical(data.choice);
PDW = logical(data.PDW);

% to avoid log(0) issues:
% minP = eps;
minP = 1e-300;
pRight_model_trialwise(pRight_model_trialwise==0) = minP; 
% pRight_model_trialwise(pRight_model_trialwise==1) = 1-minP;
pHigh_model_trialwise(pHigh_model_trialwise<=0) = minP; 
% pHigh_model_trialwise(pHigh_model_trialwise>=1) = 1-minP;

% CHOICE
% log likelihood of rightward choice on each trial, under binomial assumptions:
LL_choice = sum(log(pRight_model_trialwise(choice))) + sum(log(1-pRight_model_trialwise(~choice)));

% RT
if options.RTtask     
    
%     % likelihood of mean RTs for each coherence, under Gaussian approximation
% %     L_RT = 1./(sigmaRT_data*sqrt(2*pi)) .* exp(-(meanRT_model-meanRT_data).^2 ./ (2*sigmaRT_data.^2));
%     L_RT_high = 1./(sigmaRThigh_data*sqrt(2*pi)) .* exp(-(meanRThigh_model-meanRThigh_data).^2 ./ (2*sigmaRThigh_data.^2));
%     L_RT_low = 1./(sigmaRTlow_data*sqrt(2*pi)) .* exp(-(meanRTlow_model-meanRTlow_data).^2 ./ (2*sigmaRTlow_data.^2));
%     LL_RT = log(max(minP,L_RT_high)) + log(max(minP,L_RT_low)); % this doesn't seem right; likelihoods are greater than 1
    
    % OR assign means to trials and sum over those, to keep same order of magnitude as other LLs:
%     I = usetrs_data;
%     L_RT = 1./(sigmaRT_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RT_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRT_data_trialwise(I).^2));
%     L_RT(L_RT==0) = min(L_RT(L_RT~=0)); % this seems weird; there shouldn't be any 0 vals anyway...
%     LL_RT = nansum(log(L_RT(:))); % sum over all conditions (or trials)    
    I = usetrs_data & data.PDW_preAlpha==1;
    L_RT_high = 1./(sigmaRThigh_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RThigh_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRThigh_data_trialwise(I).^2));    
    L_RT_high(L_RT_high==0) = minP;
    I = usetrs_data & data.PDW_preAlpha==0;    
    L_RT_low = 1./(sigmaRTlow_data_trialwise(I)*sqrt(2*pi)) .* exp(-(RTlow_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRTlow_data_trialwise(I).^2));
    L_RT_low(L_RT_low==0) = minP;

    LL_RT = nansum(log(L_RT_high(:))) + nansum(log(L_RT_low(:)));

    % next: fit full distributions  

    
    
    
    
end

% CONF
if options.conftask==1 % SEP
    keyboard % this is unfinished
    % likelihood of mean conf ratings for each condition, under Gaussian approximation
    L_conf = 1./(sigmaConf_data*sqrt(2*pi)) .* exp(-(meanConf_model - meanConf_data).^2 ./ (2*sigmaConf_data.^2));
    L_conf(L_conf==0) = min(L_conf(L_conf~=0));
    LL_conf = nansum(log(L_conf(:))); % sum over all conditions
elseif options.conftask==2 % PDW
    % log likelihood of high bet on each trial, under binomial assumptions:
    LL_conf = sum(log(pHigh_model_trialwise(PDW))) + sum(log(1-pHigh_model_trialwise(~PDW)));
    
    % OR multinomial likelihood with n=1 and k=4 (categorial distribution)
    RH = data.choice==1 & data.PDW==1;
    RL = data.choice==1 & data.PDW==0;
    LH = data.choice==0 & data.PDW==1;
    LL = data.choice==0 & data.PDW==0;
    LL_multinom = sum(log(max(pRightHigh_model_trialwise(RH),minP))) + sum(log(max(pRightLow_model_trialwise(RL),minP))) + ...
                  sum(log(max(pLeftHigh_model_trialwise(LH),minP)))  + sum(log(max(pLeftLow_model_trialwise(LL),minP)));
    
%     % sanity check: go coh by coh, should be the same
%     LL_choice_pdw2 = 0;
%     for c = 1:length(cohs)
%         I = data.scoh==cohs(c);
%         RH = I & data.choice==1 & data.PDW==1;
%         RL = I & data.choice==1 & data.PDW==0;
%         LH = I & data.choice==0 & data.PDW==1;
%         LL = I & data.choice==0 & data.PDW==0;
%         LL_c = sum(log(max(pRightHigh_model_trialwise(RH),minP))) + sum(log(max(pRightLow_model_trialwise(RL),minP))) + ...
%                sum(log(max(pLeftHigh_model_trialwise(LH),minP)))  + sum(log(max(pLeftLow_model_trialwise(LL),minP)));
%         LL_choice_pdw2 = LL_choice_pdw2 + LL_c;
%         % yep, close enough.
%            
% %         % is this the same as summing over k, using the total number of
% %         % occurrences of each outcome (as in "bag of words" blog post)?
% %         % [cannot include the factorial terms, xi is too large, turns inf]
% %         LL_choice_pdw3 = ...
% %         sum(RH)*log(max(pRightHigh_model(c),minP)) + ...
% %         sum(RL)*log(max(pRightLow_model(c),minP)) + ...
% %         sum(LH)*log(max(pLeftHigh_model(c),minP)) + ...
% %         sum(LL)*log(max(pLeftLow_model(c),minP));
% %           % definitely not. oh well.
%     
%     end

end

% total -LL
% err = -(LL_choice + LL_conf + LL_RT);
% OR
err = -(LL_multinom + LL_RT);


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
