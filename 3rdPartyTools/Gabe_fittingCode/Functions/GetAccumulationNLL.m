function [NLL f] = GetAccumulationNLL(params, d, stimStdev, useChoiceFlag, genOutput, biasCorrFlag, varargin)

% Calculates the negative log-likelihood of data in "d" given parameters of
% an accumulation (drift-diffusion) model. 

% d takes the following form:

%   d = 
% 
%   1 x nStimulusConditions struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     RTs         (nTrials x 1, Measured reaction time for each trial in this stim cond.)
%     choice      (nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)
%
%
%
%   useChoiceFlag takes the value 0, 1, or 3 (don't ask)
%       if 0, calculate NLL using only RTs
%       if 1, calculate NLL using both RTs and choice data
%       if 3, calculate NLL using only choice data. Remember that you
%       cannot predict RTs with fits to choice data.



%%

tNDmax_right = 5 ;
tNDmax_left  = 5 ;

% Define params
k         = params(1) ; % Kappa
A1        = params(2) ; % Bound
tND_right = params(3) ; % non-decision time, righward
tND_left  = params(4) ; % non-decision time, leftward
offset    = params(5) ; % Bias
v         = params(6) ; 


% Hyperparameters
lRate = 0 ;
v = 0 ;
beta = 1 ;



for i=1:length(varargin)
    
    if isequal(varargin{i},'k')
        if varargin{i+1} ~= 0
            k=varargin{i+1};
        end
    end
    
    if isequal(varargin{i},'B')
        if varargin{i+1} ~= 0
            A1=varargin{i+1};
        end
    end
    
    if isequal(varargin{i},'tNDmax_right')
        if varargin{i+1} ~= 0
            tNDmax_right=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tND_right')
        if varargin{i+1} ~= 0
            tND_right=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tNDmax_left')
        if varargin{i+1} ~= 0
            tNDmax_left=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tND_left')
        if varargin{i+1} ~= 0
            tND_left=varargin{i+1};
        end
    end

    if isequal(varargin{i},'tNDmax')
        if varargin{i+1} ~= 0
            tNDmax_right=varargin{i+1};
            tNDmax_left=varargin{i+1};
        end
    end
        
end




A2 = A1 ;
tND_left = tND_right ;

% Make sure tND can't be above some maximum
if tND_right > tNDmax_right
    tND_right = tNDmax_right ;
end

if tND_left > tNDmax_left
    tND_left = tNDmax_left ;
end

% Fit glm to get bias term. Bias term will be used to select trials where
% choice matched the sign of Mu.
if biasCorrFlag

    if length([d.coherence]) > 2
        % Create predictor matrix
        for i = 1:length(d)
            temp{i}  = repmat(d(i).coherence, length(d(i).RTs), 1) ;
        end

        glmX = cat(1, temp{:}) ;                            % Vector containing the coherence on each trial
        glmY = reshape([d.choice],numel([d.choice]), 1) ;   % Vector containing responses on each trial

        % Fit glm
        Bs   = glmfit(glmX, glmY, 'binomial', 'logit') ;
        bias = Bs(1) / Bs(2) ;

    else
        bias = 0 ;
    end


    % Get reaction times for trials where choice matches Mu (i.e. rightward
    % choices when C > bias and leftward choices when C < bias). Calc mean and
    % variance.
    for i = 1:length(d) 
        if d(i).coherence >= -bias
            goodRTs{i} = d(i).RTs(d(i).choice == 1) ;
        else
            goodRTs{i} = d(i).RTs(d(i).choice == 0) ;
        end

        meanRTs(i)   = nanmean(goodRTs{i}) ;
        predRTsem(i) = sqrt(nanvar(goodRTs{i}) ./ length(goodRTs{i})) ; 
        nTrials(i) = length(goodRTs{i}) ;

    end

else
    
    meanRTs = nanmean([d.RTs]) ;
    predRTsem = sqrt(nanstd([d.RTs]) ./ [d.nTrials]) ;
    nTrials = [d.nTrials] ;
    
end
    


x = [d.coherence] - offset ;

mu   = sign(x) .* (abs(k .* x) .^ beta) ;
sigma2 = 1 + abs(x)*v ;



% Change zero mu to very small to avoid using limits
mu(mu==0) = 1e-8 ;

%%


% Calculate model predictions

    % Get accuracy prediction
%     predChoice = lRate + (1 - 2 * lRate) .*  1./(1+(exp(-2*mu*A1./sigma2))) ;
    predChoice = lRate + (1 - 2 * lRate) .*  ( (exp(2*mu*A2) - 1) ./ (exp(2*mu*A2) - exp(-2*mu*A1)) ) ; 


    % Get RT prediction
    predRTmean = zeros(1, length(mu)) ;
    predDT = zeros(1, length(mu)) ;
    leftMu = mu < 0 ;
    muL = mu(leftMu) ;
    muR = mu(~leftMu) ;


%     predDT     = (A1 ./ mu) .* tanh(mu*A1./sigma2) ; 
    predDT(leftMu)  = ( (A1 + A2) ./ muL) .* coth( (A1+A2).*muL) - (A2./muL).*coth(A2.*muL) ;
    predDT(~leftMu) = ( (A1 + A2) ./ muR) .* coth( (A1+A2).*muR) - (A1./muR).*coth(A1.*muR) ;
    predRTmean(leftMu)  = predDT(leftMu)  + tND_left ;
    predRTmean(~leftMu) = predDT(~leftMu) + tND_right ;
   



% Calculate log-likelihood of observed mean RTs, given predicted mean RTs
% and predicted sem of RTs 
likelihoodRT = normpdf(meanRTs, predRTmean, predRTsem) ;
likelihoodRT(likelihoodRT==0) = 1e-50 ;
rtNLL        = -sum(log(likelihoodRT)) ;

% Get likelihood for choice
likelihoodChoice = binopdf(nansum([d.choice]), [d.nTrials], predChoice) ;
likelihoodChoice(likelihoodChoice==0) = 1e-50 ;
choiceNLL        = -sum(log(likelihoodChoice)) ;

% Calculate total NLL
if useChoiceFlag == 1
    NLL = rtNLL + choiceNLL ;
elseif useChoiceFlag == 0
    NLL = rtNLL ;
elseif useChoiceFlag == 3
    NLL = choiceNLL ;
end

if NLL == inf || isnan(NLL)
    NLL = 1e80 ;
end

%%

if genOutput

    c4p = linspace(d(1).coherence, d(end).coherence, 500) ;
%     c4p = -1.05:.005:1.05 ;
    mu4p = (k*(c4p - offset))  ;
    v4p  = 1 + (v*(abs(c4p - offset))) ;
    

    f.NLL = NLL ;
    f.params = params ;

%     f.predChoice = lRate + (1 - 2 * lRate) .*  1./(1+(exp(-2*mu4p*A1./v4p)));
    f.predChoice = lRate + (1 - 2 * lRate) .*  ( (exp(2*mu4p*A2) - 1) ./ (exp(2*mu4p*A2) - exp(-2*mu4p*A1)) ) ; 

    f.predRTmean = zeros(1, length(mu4p)) ;
    f.predDTmean = zeros(1, length(mu4p)) ;
    lMu = mu4p < 0 ;
    mu4pL = mu4p(lMu) ;
    mu4pR = mu4p(~lMu) ;

%     f.predDTmean = ((A1 ./ mu4p) .* tanh(mu4p*A1./v4p)) ; 
    f.predDTmean(lMu)  = ( (A1 + A2) ./ mu4pL) .* coth( (A1+A2).*mu4pL) - (A2./mu4pL).*coth(A2.*mu4pL) ;
    f.predDTmean(~lMu) = ( (A1 + A2) ./ mu4pR) .* coth( (A1+A2).*mu4pR) - (A1./mu4pR).*coth(A1.*mu4pR) ;
    f.predRTmean(lMu)  = f.predDTmean(lMu)  + tND_left ;
    f.predRTmean(~lMu) = f.predDTmean(~lMu) + tND_right ;

    f.meanRTs = meanRTs ;
    f.nTrials = nTrials ;
    f.predRTsem  = predRTsem ;
    f.c4p = c4p ;
    f.sigma2 = sigma2 ;
    f.rtNLL = rtNLL ;
    
    if A1 ~= A2 || tND_right ~= tND_left
        [~, i] = closest_to(f.c4p, offset);
        f.predRTmean(i) = nan;
    end

else
    
    f = [] ;
    
end



end