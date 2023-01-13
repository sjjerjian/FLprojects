function [err,fit,parsedFit] = errfcn_DDM_2D_wConf_noMC(param, guess, fixed, data, options)

% updated model calculations to Steven's method: 
% SJ 10-11-2021 no Monte Carlo simulation for fitting, it's redundant!
% just use model predictions directly

% uses UNSIGNED cohs for MOI calculation

tic

global call_num

% retrieve the full parameter set given the adjustable and fixed parameters 
param = getParam(param, guess, fixed);
if isfield(data,'dur')
    max_dur = max(data.dur)/1000; % max stimulus duration (s)
else
    max_dur = 3;
end

k = param(1);
B = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3)); % or negative thetas
alpha = param(4); % base rate of low-conf choices
% % Tnd = param(5); % fixed Tnd, for now
if length(param)==5
    TndR = param(5); 
    TndL = param(5);
else
    TndR = param(5); % optionally, separate Tnds for L and R
    TndL = param(6);
end


% use method of images to calculate PDFs of DV and mapping to log odds corr
R.t = 0.001:0.001:max_dur;
R.Bup = B;
R.drift = k * unique(data.coherence); % takes only unsigned drift rates
R.drift_freq = hist(data.coherence,unique(data.coherence))'/length(data.coherence); % relative frequency of each coh(drift)
R.lose_flag = 1; % we always need the losing densities
R.plotflag = options.plot; % 1 = plot, 2 = plot nicer and export_fig (eg for talk)
%     R.plotflag = 0; % manual override
P = images_dtb_2d(R); % /WolpertMOI

% for fitting RT dists, need to  clip any data.RTs that may be past the
% limit in the mapping
dataRT_forPDF = data.RT;
if sum(dataRT_forPDF>P.t(end))/length(dataRT_forPDF)>0.05
    warning('Lots of long RTs, check the data or extend the mapping')
end
dataRT_forPDF(dataRT_forPDF>P.t(end)) = P.t(end);

%% calculate likelihood of the observed data given the model parameters

% usetrs_data = data.correct | data.coherence<1e-6; % use only correct (OR ~0 hdg) trials for RT fits
usetrs_data = true(size(data.correct)); % try all trials

% coh-wise vars (for parsedFit, and some likelihood calculations)
cohs = unique(data.scoh);
n = nan(length(cohs),1);
pRight_model = n;
pHigh_model = n;
pRight_High_model = n;
pRight_Low_model = n;
pHigh_Corr_model = n;
pHigh_Err_model = n;
n_RT_data = n;
n_RThigh_data = n;
n_RTlow_data = n;
meanRT_model = n;
meanRThigh_model = n;
meanRTlow_model = n;
meanRT_data = n;
meanRThigh_data = n;
meanRTlow_data = n;
% % % sigmaRT_data = n;  % obsolete
% % % sigmaRThigh_data = n;
% % % sigmaRTlow_data = n;
sigmaRT_model = n;
sigmaRThigh_model = n;
sigmaRTlow_model = n;
meanConf_model = n;
meanConf_data = n;
sigmaConf_data = n;
nCor = n;

% trialwise vars
m = nan(length(data.coherence),1);
% marginals
pRight_model_trialwise = m;
pHigh_model_trialwise  = m;
% intersections
pRightHigh_model_trialwise = m;
pRightLow_model_trialwise = m;
pLeftHigh_model_trialwise = m;
pLeftLow_model_trialwise = m;
% conditionals
pRight_High_model_trialwise = m;
pRight_Low_model_trialwise = m;
pHigh_Corr_model_trialwise = m;
pHigh_Err_model_trialwise = m;
% other
conf_model_trialwise = m;
RT_model_trialwise = m;
RThigh_model_trialwise = m;
RTlow_model_trialwise = m;
% % % sigmaRT_data_trialwise = m;  % obsolete
% % % sigmaRThigh_data_trialwise = m;
% % % sigmaRTlow_data_trialwise = m;
sigmaRT_model_trialwise = m;
sigmaRThigh_model_trialwise = m;
sigmaRTlow_model_trialwise = m;
L_RT_fromPDF = m;
L_RT_high_fromPDF = m;
L_RT_low_fromPDF = m;

for c = 1:length(cohs) % loop through signed cohs, because that is what defines a trial in the data
    
    % index for images_dtb (unsigned coh/drift corresponding to this signed coh)                    
    uc = abs(cohs(c))==(cohs(cohs>=0));  % *** marking differences in signed vs. unsigned ver ***
      
    % grab the distributions from images_dtb:
    Pxt1 = squeeze(P.up.distr_loser(uc,:,:))'; % losing race pdf given correct
    Pxt2 = squeeze(P.lo.distr_loser(uc,:,:))'; % losing race pdf given incorrect

    % ************
    % CHOICE & PDW
    % ************
    if cohs(c)<0 % leftward motion
        % if coh is negative, then pRight is p(incorrect), aka P.lo
        pRight = P.lo.p(uc)/(P.up.p(uc)+P.lo.p(uc)); % normalized by total absorbed probability 
        pLeft = 1-pRight; % (evidently, the unabsorbed probability can be ignored when it comes to pRight)

        % calculate probabilities of the four outcomes: right/left x high/low
        % these are the intersections, e.g. P(R n H) [n is used as the upside-down U symbol]
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
    Ptot = pRightHigh + pRightLow + pLeftHigh + pLeftLow;
    if abs(Ptot-1)>0.01; warning('Ptot'); fprintf('c=%d, Ptot=%f\n',c,Ptot); end
    pRightHigh = pRightHigh/Ptot;
    pRightLow = pRightLow/Ptot;
    pLeftHigh = pLeftHigh/Ptot;
    pLeftLow = pLeftLow/Ptot;

% % %     % copy intersections to trials, for multinomial likelihood
% % %     % does NOT take into account alpha!!!
% % %     Jdata = data.scoh==cohs(c); % select trials of this coh
% % %     pRightHigh_model_trialwise(Jdata) = pRightHigh;
% % %     pRightLow_model_trialwise(Jdata) = pRightLow;
% % %     pLeftHigh_model_trialwise(Jdata) = pLeftHigh;
% % %     pLeftLow_model_trialwise(Jdata) = pLeftLow;
    % this step now occurs below, after incorporating alpha 

    % by definition of conditional probability: P(A|B) = P(A n B) / P(B)
    pHigh_Right = pRightHigh / pRight;
    pLow_Right = pRightLow / pRight; 
    pHigh_Left = pLeftHigh / pLeft;
    pLow_Left = pLeftLow / pLeft;  
    
    % by law of total probability
    pHigh = pHigh_Right*pRight + pHigh_Left*pLeft;
    pLow = pLow_Right*pRight + pLow_Left*pLeft;

    % by Bayes' rule:
    pRight_High = pHigh_Right * pRight / pHigh;
    pRight_Low = pLow_Right * pRight / pLow;
    pLeft_High = pHigh_Left * pLeft / pHigh;
    pLeft_Low = pLow_Left * pLeft / pLow;

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
    
    % and use (rearranged) definition of conditional probability to adjust
    % the intersections
    pRightHigh_wAlpha = pHigh_wAlpha * pRight_High;
    pRightLow_wAlpha = (1-pHigh_wAlpha) * pRight_Low;
    pLeftHigh_wAlpha = pHigh_wAlpha * pLeft_High;
    pLeftLow_wAlpha = (1-pHigh_wAlpha) * pLeft_Low;

    % still should sum to 1
    Ptot = pRightHigh_wAlpha + pRightLow_wAlpha + pLeftHigh_wAlpha + pLeftLow_wAlpha;
    if abs(Ptot-1)>1e-6; warning('Ptot'); fprintf('c=%d, Ptot=%f\n',c,Ptot); end
    pRightHigh_wAlpha = pRightHigh_wAlpha/Ptot;
    pRightLow_wAlpha = pRightLow_wAlpha/Ptot;
    pLeftHigh_wAlpha = pLeftHigh_wAlpha/Ptot;
    pLeftLow_wAlpha = pLeftLow_wAlpha/Ptot;

    % copy intersections to trials, for multinomial likelihood
    % *NOW INCLUDES effect of alpha
    Jdata = data.scoh==cohs(c); % select trials of this coh
    pRightHigh_model_trialwise(Jdata) = pRightHigh_wAlpha;
    pRightLow_model_trialwise(Jdata) = pRightLow_wAlpha;
    pLeftHigh_model_trialwise(Jdata) = pLeftHigh_wAlpha;
    pLeftLow_model_trialwise(Jdata) = pLeftLow_wAlpha;
    
    % copy to vectors for parsedFit struct
    pRight_model(c) = pRight;
    pHigh_model(c) = pHigh_wAlpha; % includes the effect of alpha
    pRight_High_model(c) = pRight_High; % does NOT include the effect of alpha
    pRight_Low_model(c) = pRight_Low; % does NOT include the effect of alpha
    
    pHigh_Corr_model(c) = pHigh_Corr; % includes the effect of alpha
    pHigh_Err_model(c) = pHigh_Err; % includes the effect of alpha

    % copy to trials for likelihood
    pRight_model_trialwise(Jdata) = pRight_model(c);
    pHigh_model_trialwise(Jdata) = pHigh_model(c);
        % note these are *conditionals*, not intersections;
    pRight_High_model_trialwise(Jdata) = pRight_High_model(c);
    pRight_Low_model_trialwise(Jdata) = pRight_Low_model(c);
    pHigh_Corr_model_trialwise(Jdata) = pHigh_Corr_model(c);
    pHigh_Err_model_trialwise(Jdata) = pHigh_Err_model(c);
        
    % ************
    % RT
    % ************
    nCor(c) = sum(Jdata & (data.correct | data.coherence<1e-6));
    if options.RTtask
        
        % different Tnds for different choices, NOT for different cohs,
        % so we simply take a weighted average based on P(right)
        Tnd = TndR*pRight + TndL*(1-pRight);
        
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
        
        meanRT_model(c) = pHB*P.up.mean_t(uc) + (1-pHB)*max_dur + Tnd; % weighted average of model mean and max_dur (RT when unabsorbed)
 
        % for RT conditioned on wager, it seems we don't need to do any
        % scaling by Ptb; the correct distribution is captured by the
        % conditional Pxts, we just sum them, dot-prod with t and normalize
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

        % TEMP(?): compute the variances of these PDFs for the
        % fitting-by-mean RT likelihood calculations below. This may become
        % obsolete if we can use the PDFs themselves
        
        % variance = sum((x - mean).^2 .* y) 
        % where x is the domain of the PDF, y is the PDF values over the
        % domain x, and mean is the mean of the PDF, which you can
        % calculate using the formula you provided: mean = sum(y.*x). -- thanks ChatGPT!

% % %         % SJ 12-2022 Is diag is what we want here? Var(X) = E((X-E(X))^2)
% % %         P.up.var_t = diag(P.up.pdf_t * ((R.t - P.up.mean_t).^2)' ./ P.up.p);
% % %         P.lo.var_t = diag(P.lo.pdf_t * ((R.t - P.lo.mean_t).^2)' ./ P.lo.p);

                                                                % ___note the PDFs are normalized___
        sigmaRT_model(c) = sqrt(sum((R.t - P.up.mean_t(uc)).^2 .* P.up.pdf_t(uc,:)/sum(P.up.pdf_t(uc,:))));
        sigmaRThigh_model(c) = sqrt(sum((R.t - meanRThigh).^2 .* PxtAboveTheta/sum(PxtAboveTheta)));
        sigmaRTlow_model(c) = sqrt(sum((R.t - meanRTlow).^2 .* PxtBelowTheta/sum(PxtBelowTheta)));
            % (note these do not include any variance added from convolution w Tnd dist)
        sigmaRT_model_trialwise(Jdata) = sigmaRT_model(c);
        sigmaRThigh_model_trialwise(Jdata) = sigmaRThigh_model(c);
        sigmaRTlow_model_trialwise(Jdata) = sigmaRTlow_model(c);
        
        % copy to trials
        RT_model_trialwise(Jdata) = meanRT_model(c);
        RThigh_model_trialwise(Jdata) = meanRThigh_model(c);
        RTlow_model_trialwise(Jdata) = meanRTlow_model(c);
        
        
        % grab the mean and sd(se?) for corresponding data
        I = Jdata & usetrs_data;
        n_RT_data(c) = sum(I);
        meanRT_data(c) = mean(data.RT(I));
% % %         sigmaRT_data(c) = std(dataRT(I))/sqrt(sum(I)); % SD or SEM? --*obsolete even if not using full dists, use Sigma from the dist itself
% % %         sigmaRT_data_trialwise(Jdata) = sigmaRT_data(c); % copy to trials, for alternate LL calculation below
                                                               % (should this be Jdata or I?)            
            % but we can just get the likelihood directly from the PDF!
            PDF = P.up.pdf_t(uc,:)/sum(P.up.pdf_t(uc,:)); % renormalize
            TndDist = normpdf(P.t,Tnd,2e-4)/sum(normpdf(P.t,Tnd,2e-4)); 
            PDFwTnd = conv(PDF,TndDist); % convolve with ^Tnd dist (currently zero variance, approx something tiny)
            PDFwTnd = PDFwTnd(1:length(P.t))/sum(PDFwTnd(1:length(P.t))); % normalize again and trim the extra length from conv
            % figure;plot(PDF,'b');hold on; plot(PDFwTnd,'r'); % temp, just to see that the conv makes sense
            L_RT_fromPDF(I) = PDFwTnd(round(dataRT_forPDF(I)*1000))';

        % repeat for high/low
        I = Jdata & usetrs_data & data.PDW_preAlpha==1;
        n_RThigh_data(c) = sum(I);
        meanRThigh_data(c) = mean(data.RT(I));
% % %         sigmaRThigh_data(c) = std(dataRT(I))/sqrt(sum(I)); % SD or SEM? --*obsolete even if not using full dists, use Sigma from the dist itself
% % %         sigmaRThigh_data_trialwise(Jdata) = sigmaRThigh_data(c); % copy to trials, for alternate LL calculation below
            % but we can just get the likelihood directly from the PDF!
            PDF = PxtAboveTheta/sum(PxtAboveTheta);
            PDFwTnd = conv(PDF,TndDist);
            PDFwTnd = PDFwTnd(1:length(P.t))/sum(PDFwTnd(1:length(P.t)));
            L_RT_high_fromPDF(I) = PDFwTnd(round(dataRT_forPDF(I)*1000))';
            if any(isnan(L_RT_high_fromPDF(I))); keyboard; end % TEMP

        I = Jdata & usetrs_data & data.PDW_preAlpha==0;
        n_RTlow_data(c) = sum(I);
        meanRTlow_data(c) = mean(data.RT(I));
% % %         sigmaRTlow_data(c) = std(data.RT(I))/sqrt(sum(I)); % SD or SEM? --*obsolete even if not using full dists, use Sigma from the dist itself
% % %         sigmaRTlow_data_trialwise(Jdata) = sigmaRTlow_data(c); % copy to trials, for alternate LL calculation below
            % but we can just get the likelihood directly from the PDF!
            PDF = PxtBelowTheta/sum(PxtBelowTheta); 
            PDFwTnd = conv(PDF,TndDist);
            PDFwTnd = PDFwTnd(1:length(P.t))/sum(PDFwTnd(1:length(P.t)));
            L_RT_low_fromPDF(I) = PDFwTnd(round(dataRT_forPDF(I)*1000))';
            if any(isnan(L_RT_low_fromPDF(I))); keyboard; end % TEMP
            
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
    parsedFit.pRightHigh = pRight_High_model;
    parsedFit.pRightLow = pRight_Low_model;
    parsedFit.pHighCorr = pHigh_Corr_model;
    parsedFit.pHighErr = pHigh_Err_model;
end
if options.RTtask
    parsedFit.RTmean = meanRT_model;
    if options.conftask==2 % PDW
        parsedFit.RTmeanHigh = meanRThigh_model;
        parsedFit.RTmeanLow = meanRTlow_model;
    end
end


%% Next, calculate error (negative log-likelihood)

% convert data vars to logicals
choice = logical(data.choice);
PDW = logical(data.PDW);

% use this to avoid log(0) issues
minP = 1e-300;

% CHOICE
% log likelihood of rightward choice on each trial, under binomial assumptions:
pRight_model_trialwise(pRight_model_trialwise==0) = minP; 
LL_choice = sum(log(pRight_model_trialwise(choice))) + sum(log(1-pRight_model_trialwise(~choice)));

% % steven method, based on Palmer: choice proportions, not trials
% parsedData = Dots_parseData(data,options.conftask,options.RTtask,0);
% L_choice = binopdf(round(parsedData.pRight.*parsedData.n_all),parsedData.n_all,pRight_model); 
% L_choice(L_choice==0) = minP;
% LL_choice = nansum(log(L_choice(:)));

% CHOICE, from conditionals
pRight_Low_model_trialwise(pRight_Low_model_trialwise==0) = minP;
pp = pRight_Low_model_trialwise(data.PDW==0); % conditionalize first
cc = choice(data.PDW==0);
LL_choice_low = sum(log(pp(cc))) + sum(log(1-pp(~cc)));

pRight_High_model_trialwise(pRight_High_model_trialwise==0) = minP;
pp = pRight_High_model_trialwise(data.PDW==1); % conditionalize first
cc = choice(data.PDW==1);
LL_choice_high = sum(log(pp(cc))) + sum(log(1-pp(~cc)));

LL_choice_cond = LL_choice_low + LL_choice_high; % currently unused

% RT
if options.RTtask     
    
%     sigmaRT_data = sigmaRT_data^2 % ***this seems more comparable to the theoretical SEM from Palmer, keep it as an option.

%     L_RToption = 'means_nosep';
%     L_RToption = 'trials_nosep';
%     L_RToption = 'means_sep';
%     L_RToption = 'trials_sep';
    L_RToption = 'full_dist';

    switch L_RToption
        case 'means_nosep'
            
            % Option 1a: means, not sep
            % (fit mean RTs for each coh under Gaussian approximation, NOT separated by high/low bet)
            L_RT = 1./(sigmaRT_model*sqrt(2*pi)) .* exp(-(meanRT_model-meanRT_data).^2 ./ (2*sigmaRT_model.^2)) / 1000; % div by 1000 to account for units (sec/ms)
                    % *** NOT SURE ABOUT THE FACTOR OF 1000 ***
            L_RT(L_RT==0) = minP; % kluge to remove zeros, which are matlab's response to  exp(-[too large a number]), otherwise the log will give you infinity
            LL_RT = sum(log(L_RT(:))); % sum over all conditions (or trials)
            LL_RT_high = NaN; LL_RT_low = NaN;
            
        case 'trials_nosep'
    
            % Option 1b: trials, not sep
            % (assign means to trials and sum over those, to keep same order of mag as other LLs)
            I = usetrs_data;
            L_RT = 1./(sigmaRT_model_trialwise(I)*sqrt(2*pi)) .* exp(-(RT_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRT_model_trialwise(I).^2)) / 1000;
            L_RT(L_RT==0) = minP;
            LL_RT = sum(log(L_RT(:)));
            LL_RT_high = NaN; LL_RT_low = NaN;
            
        case 'means_sep'
    
            % Option 2a: means, sep hi/lo
            % (fit mean RT for each coherence, separately for high/low bet)
            L_RT_high = 1./(sigmaRThigh_model*sqrt(2*pi)) .* exp(-(meanRThigh_model-meanRThigh_data).^2 ./ (2*sigmaRThigh_model.^2)) / 1000;
            L_RT_high(L_RT_high==0) = minP;
                % WEIGHT LL by nTrials, e.g. since there are fewer low bets
                % at high coh, a miss there should penalize less
            weightedLL = n_RThigh_data/sum(n_RThigh_data) .* log(L_RT_high) * length(cohs);
            LL_RT_high = sum(weightedLL);

            L_RT_low = 1./(sigmaRTlow_model*sqrt(2*pi)) .* exp(-(meanRTlow_model-meanRTlow_data).^2 ./ (2*sigmaRTlow_model.^2)) / 1000;
            L_RT_low(L_RT_low==0) = minP;
            weightedLL = n_RTlow_data/sum(n_RTlow_data) .* log(L_RT_low) * length(cohs);
            LL_RT_low = sum(weightedLL);

            LL_RT = LL_RT_high + LL_RT_low;
            
        case 'trials_sep'
    
            % Option 2b: trials, sep hi/lo
            % (separate by high/low bet and assign to trials, to keep order of mag)
        %     I = usetrs_data & data.PDW_preAlpha==1; %??
            I = usetrs_data & data.PDW==1;
            L_RT_high = 1./(sigmaRThigh_model_trialwise(I)*sqrt(2*pi)) .* exp(-(RThigh_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRThigh_model_trialwise(I).^2)) / 1000;    
            L_RT_high(L_RT_high==0) = minP;
            LL_RT_high = sum(log(L_RT_high));
            
        %     I = usetrs_data & data.PDW_preAlpha==0; %??
            I = usetrs_data & data.PDW==0;    
            L_RT_low = 1./(sigmaRTlow_model_trialwise(I)*sqrt(2*pi)) .* exp(-(RTlow_model_trialwise(I) - data.RT(I)).^2 ./ (2*sigmaRTlow_model_trialwise(I).^2)) / 1000;
            L_RT_low(L_RT_low==0) = minP;
            LL_RT_low = sum(log(L_RT_low));
            
            LL_RT = LL_RT_high + LL_RT_low;
            
        case 'full_dist'
            
            % Option 3: fit full RT distributions            
%             L_RT_fromPDF(L_RT_fromPDF==0) = minP;
%             LL_RT = sum(log(L_RT_fromPDF)); % no need to use uncond

            if sum(~isnan(L_RT_high_fromPDF))+sum(~isnan(L_RT_low_fromPDF)) ~= length(data.RT)
                keyboard
            end
            L_RT_high_fromPDF(L_RT_high_fromPDF==0) = minP;
            LL_RT_high = nansum(log(L_RT_high_fromPDF));
            
            L_RT_low_fromPDF(L_RT_low_fromPDF==0) = minP;
            LL_RT_low = nansum(log(L_RT_low_fromPDF));
                % these match the trialwise-by-means pretty well! reassuring.
                
            LL_RT = LL_RT_high + LL_RT_low;
    end

end

% CONF
if options.conftask==1 % SEP
    keyboard % this is unfinished
    % likelihood of mean conf ratings for each condition, under Gaussian approximation
    L_conf = 1./(sigmaConf_data*sqrt(2*pi)) .* exp(-(meanConf_model - meanConf_data).^2 ./ (2*sigmaConf_data.^2));
    L_conf(L_conf==0) = min(L_conf(L_conf~=0));
    LL_conf = sum(log(L_conf(:))); % sum over all conditions
elseif options.conftask==2 % PDW
    % log likelihood of high bet on each trial, under binomial assumptions:
    LL_conf = sum(log(pHigh_model_trialwise(PDW))) + sum(log(1-pHigh_model_trialwise(~PDW)));
    
    % CONF, from conditionals
    pHigh_Corr_model_trialwise(pHigh_Corr_model_trialwise==0) = minP;
    pp = pHigh_Corr_model_trialwise(data.correct==1); % conditionalize first
    cc = PDW(data.correct==1);
    LL_conf_corr = sum(log(pp(cc))) + sum(log(1-pp(~cc)));
    
    pHigh_Err_model_trialwise(pHigh_Err_model_trialwise==0) = minP;
    pp = pHigh_Err_model_trialwise(data.correct==0); % conditionalize first
    cc = PDW(data.correct==0);
    LL_conf_err = sum(log(pp(cc))) + sum(log(1-pp(~cc)));
    
    LL_conf_cond = LL_conf_corr + LL_conf_err; % currently unused
    
    % OR multinomial likelihood with n=1 and k=4 (categorial distribution)
    RH = data.choice==1 & data.PDW==1;
    RL = data.choice==1 & data.PDW==0;
    LH = data.choice==0 & data.PDW==1;
    LL = data.choice==0 & data.PDW==0;
    % this uses the intersections, not conditionals
    LL_multinom = sum(log(max(pRightHigh_model_trialwise(RH),minP))) + sum(log(max(pRightLow_model_trialwise(RL),minP))) + ...
                  sum(log(max(pLeftHigh_model_trialwise(LH),minP)))  + sum(log(max(pLeftLow_model_trialwise(LL),minP)));

% %     % double check that they all sum to 1:
% %     pRightHigh_model_trialwise + pRightLow_model_trialwise + pLeftHigh_model_trialwise + pLeftLow_model_trialwise

end


% total -LL
% err = -(LL_choice_cond + LL_conf_cond + LL_RT);
% OR
% err = -(LL_choice_cond + LL_conf_cond);
% OR
err = -(LL_multinom + LL_RT);
     %** NOW LL_RT is greater by factor of ~6. let's keep it that way and see if it matters
% OR
% err = -LL_multinom;


% OR: fit choice and RT first, then hold those fixed to find theta
% (Miguel's method)
% err = -(LL_choice + LL_RT);


%% print progress report!
if options.feedback
    fprintf('\n\n\n****************************************\n');
    fprintf('run %d\n', call_num);
    fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\talpha= %g\n\tTndR= %g\n\tTndL= %g\n', k, B, theta, alpha, TndR, TndL);
    fprintf('err: %f\n', err);
end
if options.feedback==2 && strcmp(options.fitMethod,'fms')
    figure(options.fh);
    plot(call_num, err, '.','MarkerSize',14);
    drawnow;
end
call_num = call_num + 1;
toc


%************
% temp, diagnosing LL issues
fprintf('\nLL_choice= %g\nLL_choice_low= %g\nLL_choice_high= %g\nLL_conf= %g\nLL_conf_corr= %g\nLL_conf_err= %g\nLL_multinom= %g\nLL_RT= %g\nLL_RT_low= %g\nLL_RT_high= %g\n', ...
           LL_choice,     LL_choice_low,     LL_choice_high,     LL_conf,     LL_conf_corr,     LL_conf_err,     LL_multinom,     LL_RT,     LL_RT_low,     LL_RT_high);
%************


end


% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam ( param1 , guess , fixed )
  param2(fixed==0) = param1(fixed==0);  %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess
end
