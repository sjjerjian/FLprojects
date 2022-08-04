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
theta = abs(param(3)); %or negative thetas
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

% trialwise vars
pRight_model_trialwise = nan(length(data.coherence),1);
pHigh_model_trialwise  = nan(length(data.coherence),1);
pRightHigh_model_trialwise  = nan(length(data.coherence),1);
pRightLow_model_trialwise  = nan(length(data.coherence),1);
conf_model_trialwise = nan(length(data.coherence),1);
sigmaConf_data_trialwise = nan(length(data.coherence),1);
RT_model_trialwise = nan(length(data.coherence),1);
sigmaRT_data_trialwise = nan(length(data.coherence),1);

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
meanRT_data = n;
sigmaRT_data = n;
meanConf_model = n;
meanConf_data = n;
sigmaConf_data = n;
nCor = n;

for c = 1:length(cohs) % loop through signed cohs, because that is what defines a trial in the data
    
    % index for images_dtb (unsigned coh/drift corresponding to this signed coh)                    
    uc = abs(cohs(c))==(cohs(cohs>=0));  % *** marking differences in signed vs. unsigned ver ***
      
    % grab the distributions from images_dtb:
    Pxt1 = squeeze(P.up.distr_loser(uc,:,:))'; % losing race pdf given correct
    Pxt2 = squeeze(P.lo.distr_loser(uc,:,:))'; % losing race pdf given incorrect

    % CHOICE
    if cohs(c)<0 % leftward motion
        % if coh is leftward, then pRight is p(incorrect), aka P.lo
        pRight = P.lo.p(uc)/(P.up.p(uc)+P.lo.p(uc));
        pLeft = 1-pRight;

        % by the definition of conditional probability, P(A|B) = P(A n B) / P(B)
        % (evidently, the Pxts give us the intersections, not conditionals)
        pHigh_Right = sum(sum(Pxt2.*(P.logOddsCorrMap>=theta))) / pRight;
        pHigh_Left = sum(sum(Pxt1.*(P.logOddsCorrMap>=theta))) / pLeft;
        pLow_Right = sum(sum(Pxt2.*(P.logOddsCorrMap<theta))) / pRight; 
        pLow_Left = sum(sum(Pxt1.*(P.logOddsCorrMap<theta))) / pLeft;
        
    else % rightward motion
        % if coh is leftward, then pRight is p(correct), aka P.up
        pRight = P.up.p(uc)/(P.up.p(uc)+P.lo.p(uc));
        pLeft = 1-pRight;

        % by the definition of conditional probability, P(A|B) = P(A n B) / P(B)
        % (evidently, the Pxts give us the intersections, not conditionals)
        pHigh_Right = sum(sum(Pxt1.*(P.logOddsCorrMap>=theta))) / pRight;
        pHigh_Left = sum(sum(Pxt2.*(P.logOddsCorrMap>=theta))) / pLeft;
        pLow_Right = sum(sum(Pxt1.*(P.logOddsCorrMap<theta))) / pRight; 
        pLow_Left = sum(sum(Pxt2.*(P.logOddsCorrMap<theta))) / pLeft;
    end
    
    % by law of total probability
    pHigh = pHigh_Right*pRight + pHigh_Left*pLeft;
    pLow = pLow_Right*pRight + pLow_Left*pLeft;

    % by Bayes' rule:
    pRight_High = pHigh_Right * pRight / pHigh;
    pRight_Low = pLow_Right * pRight / pLow;
    pLeft_High = pHigh_Left * pLeft / pHigh;
% % %     pLeft_Low = pLow_Left * pLeft / pLow; % unused

    % adjust the probabilities for the base rate of low-conf bets:
    % the idea is that Phigh and Plow each get adjusted down/up in
    % proportion to how close they are to 1 or 0, respectively
    pHigh_wAlpha = pHigh - alpha*pHigh;

    % Lastly, the PDW split; use Bayes to adjust P(H|R) and P(H|L)
    pHigh_Right_wAlpha = pRight_High * pHigh_wAlpha / pRight;
    pHigh_Left_wAlpha = pLeft_High * pHigh_wAlpha / pLeft;
    if cohs(c)<0 % for leftward,
        pHigh_Corr = pHigh_Left_wAlpha;
        pHigh_Err = pHigh_Right_wAlpha;
    else % for rightward,
        pHigh_Corr = pHigh_Right_wAlpha;
        pHigh_Err = pHigh_Left_wAlpha;    
    end
    
    % copy to vectors for parsedFit
    pRight_model(c) = pRight;
    pHigh_model(c) = pHigh_wAlpha; % this is the only one that includes effect of alpha
    pRightHigh_model(c) = pRight_High;
    pRightLow_model(c) = pRight_Low;
    pHighCorr_model(c) = pHigh_Corr;
    pHighErr_model(c) = pHigh_Err;

    % copy to trials
    Jdata = data.scoh==cohs(c); % select trials of this coh
    pRight_model_trialwise(Jdata) = pRight_model(c);
    pHigh_model_trialwise(Jdata) = pHigh_model(c);
    pRightHigh_model_trialwise(Jdata) = pRightHigh_model(c);
    pRightLow_model_trialwise(Jdata) = pRightLow_model(c);
    
    
    % RT
    nCor(c) = sum(Jdata & (data.correct | data.coherence<1e-6));
    if options.RTtask            
        % model RT needs to be adjusted to account for unabsorbed
        % probability (trials that don't hit the bound, usually only at low
        % cohs, when bound is high, etc). Turns out the calculation of
        % P.up.mean_t only applies to absorbed prob (bound-crossings), so
        % we simply adjust it by averaging mean_t with maxdur, weighted by
        % one minus the probability of bound crossing
        
        pHB = P.up.p(uc)+P.lo.p(uc); % probability hit bound
        meanRT_model(c) = pHB.*P.up.mean_t(uc) + (1-pHB)*max_dur + Tnd; % weighted average
        
        
        % what about RT conditioned on wager? that's more of a challenge.
        
%         % first lay out some basics. meanRT is dot product of its pdf and
%         % the time axis vector. from images_dtb:
%         P.up.mean_t = P.up.pdf_t * R.t'  ./ P.up.p;
%         P.lo.mean_t = P.lo.pdf_t * R.t'  ./ P.lo.p;

        % in the 1D code (which works!), mean RT given bound crossing AND
        % high bet (the conditional, not intersection) is:
%       RThit_high = sum((Ptb(tHi,1)+Ptb(tHi,2)) .* t_ms(tHi)) / Phit_high;

%       In that equation, tHi is the time range at which bound crossing is
%       guaranteed to result in a high bet, i.e. any times earlier than
%       where the theta criterion crosses the bound. In the 2D model, bound
%       crossing is never a guarantee but is always associated with some
%       probability of losing race being above theta; seems like the
%       analogy is useful.

%       In words, the 1D equation to me looks like the conditional
%       expectation of RT (given bound crossing, which almost always
%       happens; we can worry about weighted average with max_dur later),
%       divided by the intersection of P(BoundHit|high).
        
%       for the numerator (conditional expecation) we need conditional PDF,
        % which is: P(RT|H) = P(RT,H) / P(H)
        
        % so now we need the joint PDF P(RT,H)
        % we know how to calculate the intersection P(Right,H):
        % it's sum(sum(Pxt>theta))
        % for the joint PDF, it should be:
        % P(RT=r|H) * P(H) OR P(H|RT=r) * P(RT=r)
        
        % all of this is just moving things around! writing things in terms
        % of other things I cannot calculate. There must be an expression
        % for the P(RT=r|H) or the intersection P(RT n H), using the
        % ingredients already calculated in images_dtb (what I'm calling
        % Pxt and Ptb).

        % the rest of this is just scratch-pad...
%         Ptb1 = P.up.pdf_t(uc,:);
%         Ptb2 = P.lo.pdf_t(uc,:);
%         Pxt1Hi = sum(Pxt1.*(P.logOddsCorrMap>=theta))
%         sum(sum(Ptb1)*Pxt1Hi .* R.t) / pHigh; % nope, incorrect
                
        % condExp = ???;
        % intersect = sum(sum(Pxt1.*(P.logOddsCorrMap>=theta)));        
%         meanRThigh_model(c) = condExp/intersect;

%         RThit_high = sum((Ptb1+Ptb2) .* R.t) / (pHigh*2);
%         meanRThigh_model(c) = pHB*RThit_high + (1-pHB)*max_dur + Tnd; % NOPE!

        RT_model_trialwise(Jdata) = meanRT_model(c); % copy the mean to each trial, for 'fit' struct
        meanRT_data(c) = mean(data.RT(Jdata & usetrs_data));                 % mean and sigma for the
        sigmaRT_data(c) = std(data.RT(Jdata & usetrs_data)) / sqrt(nCor(c)); % data we want to fit
        sigmaRT_data_trialwise(Jdata) = sigmaRT_data(c); % copy the SD to each trial, for alternate LL calculation below
    end
    
end

% TEMP, placeholder until we get the correct conditionals
meanRThigh_model = meanRT_model;
meanRTlow_model = meanRT_model;

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
