function [ realChoices ] = plotChoiceTrials( choice )

removals = 0;

for l = 1:length(choice)
    if isnan(choice(l - removals)) || choice(l - removals) == 0
        % if either of those fail to report, delete the trial
        choice(l - removals) = [];
        removals = removals + 1; % removals accounts for shifts in indexes
    end
end

realChoices = choice;
% rightSelected = 0;
% rightwardProp = [];
% 
% for p = 1:length(realChoices)
%     if realChoices(p) == 2
%         rightSelected = rightSelected + 1
%     end
%     rightwardProp(p) = rightSelected/p
% end

end
