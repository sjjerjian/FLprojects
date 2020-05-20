clc
clear

% parsing data from practice dots data
% HK 10/19/18
% HK updated on 1/23 to analyze bandit data
% HK updated on 5/20 to incorporate framework Q-learning model

%% Extract data
folder = './banditData/';
subjectList = {'miguel'}; % change according to subject, case sensitive
dateRange = 20190422; % change according to date
setupFile = 'BanditTask_setup';
masterData = struct; % initialize structure for all data

for i = 1:length(subjectList) % individually pass subjects to fxn
    % copy new data to a temp struct
    tempData = Bandit_DataExtract(folder, subjectList{i}, dateRange, setupFile);
    for k = 1:length(tempData.data)
        kIndex(i) = k; % overwritten until the last run of k for i
        % copy data, conditions and cond names to a new row 'i' in masterData
        masterData.data{i,k} = tempData.data{1,k};
        masterData.conditions{i,k} = tempData.conditions{1,k};
        masterData.conditionNames{i,k} = tempData.conditionNames{1,k};
        masterData.initialParametersMerged.reward = tempData.initialParametersMerged.reward;
    end
end

%% Analyze data - general
% ok, we have our data in masterData
% let's do some things with it

[choice, accuracy, reward, coherence, highReward, lowReward] = organizeData(masterData, i, kIndex); % parse masterData
[choiceY] = plotChoiceTrials(choice);
[choiceX, rewardY, meanLeftReward, meanRightReward, stdLeftReward, stdRightReward] = plotChoiceReward(choice, reward);
[qValueLeft,qValueRight,pL,pR,alpha,beta,totalRewardSubject] = simulateRescorlaWagner(choiceX,rewardY,highReward);
[qValueLeftMA,qValueRightMA,pLMA,pRMA,totalRewardMA] = RescorlaWagnerMA(choiceX);

%% Behavioral data plotting
% Plot trials vs. choices
trialX = 1:length(choiceY);
figure(1)
subplot(2,2,1)
hold on
scatter(trialX, choiceY)
title('Behavioral trial number vs. choice')
xlabel('Trial number') 
ylabel('Choice (left vs. right)')
set(gca,'YTick',[1 2],'YTickLabels',{'A','B'})
hold off


%% Plot choices vs. reward
SELeftReward = stdLeftReward/sqrt(length(choiceX));
SERightReward = stdRightReward/sqrt(length(choiceX));
subplot(2,2,2)
hold on
plot(1, meanLeftReward, 'or')
errorbar(1, meanLeftReward, SELeftReward)
plot(2, meanRightReward, 'or')
errorbar(2, meanRightReward, SERightReward)
title('Encountered mean rewards of choices')
xlabel('Choice made') 
ylabel('Mean reward') 
set(gca,'XTick',[1 2],'XTickLabels',{'A','B'})
hold off

%% Plot q-values over trial
subplot(2,2,3)
hold on
xlength = 1:length(qValueLeft);
plot(xlength,qValueLeft,'-or','DisplayName','Left Q-value')
plot(xlength,qValueRight,'-sc','DisplayName','Right Q-value')
line([0,length(qValueLeft)],[meanLeftReward,meanLeftReward],'Color','red','LineStyle','--','DisplayName','Mean left reward')
line([0,length(qValueRight)],[meanRightReward,meanRightReward],'Color','cyan','LineStyle','--','DisplayName','Mean right reward')
title('Q-values over trial')
xlabel('Trial number') 
ylabel('Q-value')
str = {'\alpha = 0.1, \beta = 5'};
text(5,0.95,str)
legend('show')
hold off

%% Softmax probability of each choice 
subplot(2,2,4)
hold on
xlength = 1:length(pL);
plot(xlength,pL,'-or','DisplayName','Probability of left choice')
plot(xlength,pR,'-sc','DisplayName','Probability of right choice')
title('Softmax choice probabilities')
xlabel('Trial number')
ylabel('P(choice)')
str = {'\alpha = 0.1, \beta = 5'};
text(5,0.7,str)
legend('show')
hold off

%% Model-generated data plotting
% Plot q-values over trial
figure(2)
subplot(1,2,1)
hold on
xlength = 1:length(qValueLeftMA);
plot(xlength,qValueLeftMA,'-or','DisplayName','Left Q-value')
plot(xlength,qValueRightMA,'-sc','DisplayName','Right Q-value')
% line([0,length(qValueLeftMA)],[meanLeftReward,meanLeftReward],'Color','red','LineStyle','--','DisplayName','Mean left reward')
% line([0,length(qValueRightMA)],[meanRightReward,meanRightReward],'Color','cyan','LineStyle','--','DisplayName','Mean right reward')
title('Q-values over trial')
xlabel('Trial number') 
ylabel('Q-value')
str = {'\alpha = 0.1, \beta = 5'};
text(5,0.95,str)
legend('show')
hold off

% Softmax probability of each choice 
subplot(1,2,2)
hold on
xlength = 1:length(pLMA);
plot(xlength,pLMA,'-or','DisplayName','Probability of left choice')
plot(xlength,pRMA,'-sc','DisplayName','Probability of right choice')
title('Softmax choice probabilities')
xlabel('Trial number')
ylabel('P(choice)')
str = {'\alpha = 0.1, \beta = 5'};
text(5,0.7,str)
legend('show')
hold off

% Output total rewards to compare behavior vs. model
fprintf('Subject''s total reward: %1.0f \n',totalRewardSubject')
fprintf('Model''s total reward: %1.0f \n',totalRewardMA')