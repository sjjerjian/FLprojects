function [rewardGenerated] = rewardSimulate(aMA,rewardProbs,noiseFactor)

rewardGenerated = [];

if aMA == 1
    rewardGenerated = rand < (rewardProbs(1))*noiseFactor;
elseif aMA == 2
    rewardGenerated = rand < (rewardProbs(2))*noiseFactor;
end
    
end

