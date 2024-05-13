function [ choiceX, rewardY, meanLeftReward, meanRightReward, stdLeftReward, stdRightReward ] = plotChoiceReward( choice, reward )

removals = 0;

for l = 1:length(choice)
    if isnan(choice(l - removals)) || isnan(reward(l - removals)) || choice(l - removals) == 0
        % if either of those fail to report, or are a failed trial, delete the trial
        choice(l - removals) = [];
        reward(l - removals) = [];
        removals = removals + 1; % removals accounts for shifts in indexes
    end
end

rewardY = reward;
choiceX = choice;

ucho = unique(choice);
for c = 1:length(ucho) 
    meanReward(c) = mean(reward(choice==ucho(c)));
    stdReward(c) = std(reward(choice==ucho(c)));
end

meanLeftReward = meanReward(1);
meanRightReward = meanReward(2);
stdLeftReward = stdReward(1);
stdRightReward = stdReward(2);

% for k = 1:length(choice)
%     if choice(k) == 1
%         leftReward = reward(k);
%         allLeftRewards = allLeftRewards + leftReward;
%     elseif choice(k) == 2
%         rightReward = reward(k);
%         allRightRewards = allRightRewards + rightReward;
%     end
% end
% 
% meanLeftReward = mean(allLeftRewards);
% meanRightReward = mean(allRightRewards);
end