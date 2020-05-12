function [qValueLeft,qValueRight,pL,pR,totalRewardMA] = RescorlaWagnerMA(choiceX)
% based on Ten Simple Rules for Modeling by Wilson and Collins
% "simulate_M3RescorlaWagner_v1"

% simulates model data for full association.
% full association meaning any +0 reward is attributed to the underlying
% stimulus's reward rate and not attributed to an incorrect perceptual choice
% i.e. always "associated" with the stimulus

Q = [0.5 0.5]; % initial Q values - model-generated (model association)
pL = zeros(1,length(choiceX));
pR = zeros(1,length(choiceX));
qValueLeft = zeros(1,length(choiceX));
qValueRight = zeros(1,length(choiceX));

rewardProbs = [0.50 0.75]; % when right reward is higher -- rewards seem to be broken right now though
% [(Left target reward rate) (Right target reward rate)]
r = [];

% this statement isn't working for some reason -- HK 5/20
% if highReward == 1
%     rewardProbs = [0.75 0.50];
% elseif highReward == 2
%     rewardProbs = [0.50 0.75];
% else
%     fprintf('something isn''t right')
% end

% these values are subject to adjustment - but reasonable estimates for now 
alpha = 0.1;
beta = 5;
noiseFactor = 0.9; % this is a complete guess right now, will be computed from behav. data

for t = 1:length(choiceX)
    % compute choice probabilities according to softmax fxn
    p = exp(beta*Q) / sum(exp(beta*Q));
    pL(t) = p(1);
    pR(t) = p(2);
    
    % make choice according to choice probabilities
    a(t) = choose(p); % model-predicted choice - MA = model association
    
    % generate reward based on choice
    r(t) = rewardSimulate(a(t), rewardProbs, noiseFactor);
   
    % update Q-values for both behavior and model 
    delta = r(t) - Q(a(t));
    Q(a(t)) = Q(a(t)) + alpha * delta;
    qValueLeft(t) = Q(1);
    qValueRight(t) = Q(2);
    
end

totalRewardMA = sum(r);

end