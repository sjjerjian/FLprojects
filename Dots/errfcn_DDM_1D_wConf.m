function [err,fit,parsedFit] = errfcn_DDM_1D_wConf(param, guess, fixed, data, options)
    
tic

global call_num

param = getParam(param, guess, fixed);

% obsolete, I think
% % % if ~isfield(data,'PDW_preAlpha')
% % %     data.PDW_preAlpha = data.PDW;
% % % end

k = param(1);
B = abs(param(2)); % don't accept negative bound heights
theta = abs(param(3)); % no negative thetas
alpha = max([0 param(4)]); % base rate of low-conf choices (bounded at zero)
Tnd = max([0 param(5)]); % non-decision time (ms)

theta2 = theta; % allows expanded model with separate thetas for each choice (currently unused)
sigma = 1; % not a free param

%%%% scriptified this for the time being, to standardize between sim and fit
makeLogOddsMap_1D
%%%%

% accumulate log likelihood
expectedPright = nan(size(data.scoh));
expectedPrightHigh = nan(size(data.scoh));
expectedPrightLow = nan(size(data.scoh));
expectedPhigh = nan(size(data.scoh));
expectedPhighCorr = nan(size(data.scoh));
expectedPhighErr = nan(size(data.scoh));
expectedRT = nan(size(data.scoh));
expectedRThigh = nan(size(data.scoh));
expectedRTlow = nan(size(data.scoh));
t_ms = (1:size(Ptb,1))';
pRight_model = nan(length(coh_set),1);
pHigh_model = nan(length(coh_set),1);
pHighCorr_model = nan(length(coh_set),1);
pHighErr_model = nan(length(coh_set),1);
pRightHigh_model = nan(length(coh_set),1);
pRightLow_model = nan(length(coh_set),1);
meanRT_model = nan(length(coh_set),1);
meanRThigh_model = nan(length(coh_set),1);
meanRTlow_model = nan(length(coh_set),1);
try
    n = 0;
    LL = 0;
    for c = 1 : length(coh_set)

        Ptb = Ptb_coh(:,:,c); % time * bound
        Pxt = Pxt_coh(:,:,c); % xmesh * time
        
        % calculate probabilities of the four outcomes: right/left x high/low,
        % as a function of time (dur); these are the intersections, e.g. P(R n H)
        pRightHigh = cumsum(Ptb(:,2).*(bet_high_tb(:,2)==1)) + ...
                    sum(Pxt(xmesh>0,:).*(bet_high_xt(xmesh>0,:)==1),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==1))'; % half the probability of dv=0 (assuming a 50/50 guess when that happens)
        pRightLow = cumsum(Ptb(:,2).*(bet_high_tb(:,2)==0)) + ...
                    sum(Pxt(xmesh>0,:).*(bet_high_xt(xmesh>0,:)==0),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==0))';
        pLeftHigh = cumsum(Ptb(:,1).*(bet_high_tb(:,1)==1)) + ...
                    sum(Pxt(xmesh<0,:).*(bet_high_xt(xmesh<0,:)==1),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==1))';
        pLeftLow =  cumsum(Ptb(:,1).*(bet_high_tb(:,1)==0)) + ...
                    sum(Pxt(xmesh<0,:).*(bet_high_xt(xmesh<0,:)==0),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==0))';

        %ensure total prob sums to one
        Ptot = pRightHigh + pRightLow + pLeftHigh + pLeftLow;
        if any(abs(Ptot-1)>1e-3)
            warning('the probabilities do not add up to one!\n\tPtot=%f (for c= %d)\n', min(Ptot), c);
            pRightHigh = pRightHigh./Ptot;
            pRightLow = pRightLow./Ptot;
            pLeftHigh = pLeftHigh./Ptot;
            pLeftLow = pLeftLow./Ptot;
        end

        % adjust the probabilities for the base rate of low-conf bets:
        % the idea is that Phigh and Plow each get adjusted down/up in
        % proportion to how close they are to 1 or 0, respectively
        Phigh = pRightHigh + pLeftHigh;
        Phigh_wAlpha = Phigh - alpha*Phigh;
        Plow = pRightLow + pLeftLow;
        Plow_wAlpha = 1 - Phigh_wAlpha;

        % now recover the corresponding right/left proportions
        PrightHigh_wAlpha = (pRightHigh./Phigh) .* Phigh_wAlpha;
        PrightLow_wAlpha = (pRightLow./Plow) .* Plow_wAlpha;
        PleftHigh_wAlpha = (pLeftHigh./Phigh) .* Phigh_wAlpha;
        PleftLow_wAlpha = (pLeftLow./Plow) .* Plow_wAlpha;
        
        % calculate expected p(right) and P(high), given the durs
        I = data.scoh==coh_set(c);
        dur = data.dur(I);
        expectedPright(I) = (pRightHigh(dur) + pRightLow(dur)) ./ (pRightHigh(dur) + pRightLow(dur) + pLeftHigh(dur) + pLeftLow(dur));
            pRight_model(c) = mean(expectedPright(I)); % taking the mean feels wrong here.. (shouldn't matter for RT, where dur is always max_dur)
        expectedPhigh(I) = Phigh_wAlpha(dur) ./ (pRightHigh(dur) + pRightLow(dur) + pLeftHigh(dur) + pLeftLow(dur)); % only this one should take alpha into account [this is weird; see note in simDDM_postProcessing]
            pHigh_model(c) = mean(expectedPhigh(I));
        expectedPrightHigh(I) = pRightHigh(dur) ./ (pRightHigh(dur) + pLeftHigh(dur)); % conditional (P(right|high)), not intersection (PrightHigh)
            pRightHigh_model(c) = mean(expectedPrightHigh(I));
        expectedPrightLow(I) = pRightLow(dur) ./ (pRightLow(dur) + pLeftLow(dur));
            pRightLow_model(c) = mean(expectedPrightLow(I));
                        
        % pHigh conditioned on correct/error; uses the 'withAlpha' versions in the numerator
        if coh_set(c)<0 % leftward motion
            expectedPhighCorr(I) = PleftHigh_wAlpha(dur) ./ (pLeftHigh(dur) + pLeftLow(dur));
            expectedPhighErr(I) = PrightHigh_wAlpha(dur) ./ (pRightHigh(dur) + pRightLow(dur));
        else % rightward motion
            expectedPhighCorr(I) = PrightHigh_wAlpha(dur) ./ (pRightHigh(dur) + pRightLow(dur));
            expectedPhighErr(I) = PleftHigh_wAlpha(dur) ./ (pLeftHigh(dur) + pLeftLow(dur));
        end            
        pHighCorr_model(c) = mean(expectedPhighCorr(I));
        pHighErr_model(c) = mean(expectedPhighErr(I));

        if options.RTtask        
            % compute mean RT
            dI = double(coh_set(c)>=0)+1; % index for this dir (L=1, R=2)
            tHi = 1:TBint(dI)-1; % time range over which bound crossing leads to high bet
            tLo = TBint(dI):size(Ptb,1); % time range over which bound crossing leads to low bet
            % TBint is 'theta-bound intersection', see makeLogOddsMap_1D

            % compute weighted average: P(hit) * <P(RT|hit)> + P(~hit) * <P(RT|~hit)> (the latter is simply maxdur)
            if options.ignoreUnabs
                pHB = 1;
%                 pHBhigh = 1;      % this fails! why? it works in 2D...
%                 PnoHB_high = 0;   
%                 pHBlow = 1;       
%                 PnoHB_low = 0;    
            else
                pHB = sum(sum(Ptb)); % prob of crossing either bound
%                 pHBhigh = sum(Ptb(tHi,1)+Ptb(tHi,2));
%                 PnoHB_high = sum(Pxt(xmesh>0,end).*(bet_high_xt(xmesh>0,end)==1),1) + sum(Pxt(xmesh<0,end).*(bet_high_xt(xmesh<0,end)==1),1); % the portion of PRightHigh and PLeftHigh within the bounds
%                 pHBlow = sum(Ptb(tLo,1)+Ptb(tLo,2));
%                 PnoHB_low = sum(Pxt(xmesh>0,end).*(bet_high_xt(xmesh>0,end)==0),1) + sum(Pxt(xmesh<0,end).*(bet_high_xt(xmesh<0,end)==0),1); % the portion of PRightLow and PLeftLow within the bounds
            end
                        
            % calc mean RT as dot-product of PDF
            RThb = (sum((Ptb(:,1)+Ptb(:,2)) .* t_ms) / pHB);
            RTnoHB = max_dur; 
            RT_weightedAvg = pHB*RThb + (1-pHB)*RTnoHB;
            expectedRT(I) = RT_weightedAvg/1000 + Tnd; % change to s and add Tnd
            meanRT_model(c) = mean(expectedRT(I));

            % repeat for high and low bets
            pHBhigh = sum(Ptb(tHi,1)+Ptb(tHi,2));
            PnoHB_high = sum(Pxt(xmesh>0,end).*(bet_high_xt(xmesh>0,end)==1),1) + sum(Pxt(xmesh<0,end).*(bet_high_xt(xmesh<0,end)==1),1); % the portion of PRightHigh and PLeftHigh within the bounds
            RThb_high = sum((Ptb(tHi,1)+Ptb(tHi,2)) .* t_ms(tHi)) / pHBhigh;
            RT_weightedAvg_high = (pHBhigh*RThb_high + PnoHB_high*RTnoHB) / (pHBhigh + PnoHB_high); % this time, weighted avg is normalized by total P(High), which is of course not 1
            expectedRThigh(I) = RT_weightedAvg_high/1000 + Tnd; % change to s and add Tnd
            meanRThigh_model(c) = mean(expectedRThigh(I));

            pHBlow = sum(Ptb(tLo,1)+Ptb(tLo,2));
            PnoHB_low = sum(Pxt(xmesh>0,end).*(bet_high_xt(xmesh>0,end)==0),1) + sum(Pxt(xmesh<0,end).*(bet_high_xt(xmesh<0,end)==0),1); % the portion of PRightLow and PLeftLow within the bounds
            RThb_low = sum((Ptb(tLo,1)+Ptb(tLo,2)) .* t_ms(tLo)) / pHBlow;
            RT_weightedAvg_low = (pHBlow*RThb_low + PnoHB_low*RTnoHB) / (pHBlow + PnoHB_low); % this time, weighted avg is normalized by total P(Low), which is of course not 1
            expectedRTlow(I) = RT_weightedAvg_low/1000 + Tnd; % change to s and add Tnd
            meanRTlow_model(c) = mean(expectedRTlow(I));
        end

        % update log-likelihood
        dur_rightHigh = data.dur(I & data.choice==1 & data.PDW==1);
        dur_rightLow = data.dur(I & data.choice==1 & data.PDW==0);
        dur_leftHigh = data.dur(I & data.choice==0 & data.PDW==1);
        dur_leftLow = data.dur(I & data.choice==0 & data.PDW==0);
        
        % multinomial with n=1, k=4 (categorial distribution), but missing a term?
        LL_choice_pdw = sum(log(max(PrightHigh_wAlpha(dur_rightHigh),eps))) + sum(log(max(PrightLow_wAlpha(dur_rightLow),eps))) + ...
                        sum(log(max(PleftHigh_wAlpha(dur_leftHigh),eps)))   + sum(log(max(PleftLow_wAlpha(dur_leftLow),eps)));
                    
        if options.RTtask        
            % for RT, use Palmer et al. 2005, Eq. 3 & A.34+A.35 (assume zero variance in Tr aka Tnd)
            mu = abs(k*coh_set(c)); % this is mu-prime in Palmer
            if coh_set(c)==0
                VarRT = 2/3 * B^4; % Eq. limiting case for coh=0
            else
                VarRT = (B*tanh(B*mu) - B*mu.*sech(B*mu)) ./ mu.^3;
            end
            sigmaRT = sqrt(VarRT./sum(I)); % equation in text above Eq 3

            % calc separately for high bet first, 
            J = I & data.PDW_preAlpha==1;
            meanRTmodel = meanRThigh_model(c)*1000;
            meanRTdata = mean(data.RT(J))*1000;
            L_RT_high = 1./(sigmaRT*sqrt(2*pi)) .* exp(-(meanRTmodel-meanRTdata).^2 ./ (2*sigmaRT.^2));
            L_RT_high = max(L_RT_high,eps);

            % then low bet:
            J = I & data.PDW_preAlpha==0;
            meanRTmodel = meanRTlow_model(c)*1000;
            meanRTdata = mean(data.RT(J))*1000;
            L_RT_low = 1./(sigmaRT*sqrt(2*pi)) .* exp(-(meanRTmodel-meanRTdata).^2 ./ (2*sigmaRT.^2));
            L_RT_low = max(L_RT_low,eps);

            % total LL for RT
            LL_RT = log(L_RT_high)+log(L_RT_low);

            % total LL for choice, pdw, and RT for this coherence
            LL_thisC = LL_choice_pdw + LL_RT;
        else
            % total LL for choice, pdw, and RT for this coherence
            LL_thisC = LL_choice_pdw;            
        end
        
        % accumulate LL across cohs
        LL = LL + LL_thisC;
              
        % keep track of n
        n = n + length(dur_rightHigh) + length(dur_rightLow) + length(dur_leftHigh) + length(dur_leftLow);
    end
    
    % copy trial params to data struct 'fit', then replace data with model vals
    fit = data;
    fit = rmfield(fit,'choice');
    try fit = rmfield(fit,'correct'); end %#ok<TRYNC>
    fit.pRight = expectedPright;
    if options.conftask==1 % SEP
        keyboard % not ready
    %     fit.conf = conf_model_trialwise;
    %     fit.PDW = nan(size(fit.conf));
    elseif options.conftask==2 % PDW
        fit.conf = expectedPhigh;
        fit.PDW = expectedPhigh; % legacy
        fit.pRightHigh = expectedPrightHigh;
        fit.pRightLow = expectedPrightLow;
    end
    if options.RTtask            
        fit.RT = expectedRT;
        fit.RThigh = expectedRThigh;
        fit.RTlow = expectedRThigh;
    end
    
    % also stored the 'parsed' values, for later plotting
    parsedFit = struct();
    parsedFit.pRight = pRight_model;
    if options.conftask==1 % SEP
        keyboard % not ready
%         parsedFit.confMean = meanConf_model;
    elseif options.conftask==2 % PDW
        parsedFit.pHigh = pHigh_model;        
        parsedFit.pHighCorr = pHighCorr_model;
        parsedFit.pHighErr = pHighErr_model;
        parsedFit.pRightHigh = pRightHigh_model;
        parsedFit.pRightLow = pRightLow_model;
    end
    if options.RTtask
        parsedFit.RTmean = meanRT_model;
        parsedFit.RTmeanHigh = meanRThigh_model;
        parsedFit.RTmeanLow = meanRTlow_model; 
    end

catch Error
    fprintf('\n\nERROR!!!\n'); fprintf('run %d\n', call_num);
    fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\talpha= %g\nTnd= %g\n', k, B, theta, alpha, Tnd);
    fprintf('c= %d\n', c);
    for i = 1 : length(Error.stack)
        fprintf('entry %d of error stack\n', i);
        fprintf('\tfile: %s\n', Error.stack(i).file);
        fprintf('\tfunction: %s\n', Error.stack(i).name);
        fprintf('\tline: %d\n', Error.stack(i).line);
    end
    rethrow(Error);
end

err = -LL;

%print progress report
if options.feedback
    fprintf('\n\n\n****************************************\n');
    fprintf('run %d\n', call_num);
    fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\talpha= %g\n\tTnd= %g\n', k, B, theta, alpha, Tnd);
    fprintf('err: %f\n', err);
    fprintf('number of processed trials: %d\n', n);
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
function param2 = getParam(param1, guess, fixed)
  param2(fixed==0) = param1(fixed==0);       %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess
end

