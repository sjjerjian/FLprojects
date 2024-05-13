function [qValueLeft,qValueRight,pL,pR,alpha,beta,totalRewardSubject] = simulateRescorlaWagner(choiceX,rewardY,highReward)
% based on Ten Simple Rules for Modeling by Wilson and Collins
% "simulate_M3RescorlaWagner_v1"

% fits behavioral data and generates equivalent model-based parameters
% probabilities are for show and do not represent the choices made

Q = [0.5 0.5]; % initial Q values
pL = zeros(1,length(choiceX));
pR = zeros(1,length(choiceX));
qValueLeft = zeros(1,length(choiceX));
qValueRight = zeros(1,length(choiceX));

% rewardProbs = [0.50 0.75]; % in the case of highReward = 1... hardcoded for now
r = []; % reward

% these values are subject to adjustment - but reasonable estimates for now 
alpha = 0.1;
beta = 5;

for t = 1:length(choiceX)
    % compute choice probabilities according to softmax fxn
    p = exp(beta*Q) / sum(exp(beta*Q));
    pL(t) = p(1);
    pR(t) = p(2);
    
    % make choice according to choice probabilities
    a(t) = choiceX(t); % behavioral choice
    
    % generate reward based on choice
    r(t) = rewardY(t);
   
    % update Q-values for both behavior and model 
    delta = r(t) - Q(a(t));
    Q(a(t)) = Q(a(t)) + alpha * delta;
    qValueLeft(t) = Q(1);
    qValueRight(t) = Q(2);
end

totalRewardSubject = sum(r);

end
