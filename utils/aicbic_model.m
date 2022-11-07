function [aic,bic,aicc] = aicbic_model(logL,numParam,numObs)

aic = -2*logL + 2*numParam;

if nargin==3
    bic = -2*logL + log(numObs)*numParam;
    aicc = aic + (2*numParam*(numParam + 1))/(numObs - numParam - 1); % 'corrected' AIC
else
    disp('BIC/corrected AIC require numObs argument')
end