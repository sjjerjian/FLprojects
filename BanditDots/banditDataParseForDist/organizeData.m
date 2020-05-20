function [ masterChoice, masterCorrect, masterReward, masterCoherence, highReward, lowReward ] = organizeData( masterData, i, kIndex )

% fetch RTs, accuracy, coherence and the script-determined high reward
% here we lose per-subject organization, but that's changeable
masterChoice = [];
masterCorrect = [];
masterReward = [];
masterCoherence = [];

for k = 1:i % i is used here bc it still reads the height of the struct
    for j = 1:kIndex(k) % kIndex holds row lengths, see line 23 of wrapper
        try
            choice(j) = masterData.data{k,j}.behavior.choice; % these will get overwritten
        catch
            choice(j) = NaN;
        end
        try
            accuracy(j) = masterData.data{k,j}.behavior.correct;
        catch
            accuracy(j) = NaN;
        end
        try
            reward(j) = masterData.data{k,j}.reward.rewardGained;
        catch
            reward(j) = NaN;
        end
        try
            coherence(j) = masterData.conditions{k,j}.stimulus.coherence;
        catch
            coherence(j) = NaN;
        end
    end
    
    if k == 1
        masterChoice = choice;
        masterCorrect = accuracy;
        masterReward = reward;
        masterCoherence = coherence;
    else
        masterChoice = horzcat(masterChoice, choice);
        masterCorrect = horzcat(masterCorrect, accuracy);
        masterReward = horzcat(masterReward, reward);
        masterCoherence = horzcat(masterCoherence, coherence);
    end

% not sure if these values are even correct? HK 5/20
highReward = masterData.initialParametersMerged.reward.high;
lowReward = masterData.initialParametersMerged.reward.low;
    
end