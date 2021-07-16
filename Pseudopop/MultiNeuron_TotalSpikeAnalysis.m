%Miguel Vivar-Lazo
%2/27/2021
%Fetsch Lab: Mapping Firing Rates (or Total Spikes) in a Trial for
%Different Neurons
function [spikeCountRight, spikeCountLeft, selectivity_, timeMotionRight, timeMotionLeft] = MultiNeuron_TotalSpikeAnalysis(data, presentChoice, presentPDW, cohMarker, doWeNorm, wantFiringRates, individualBlocks, plotFlag)
%% Parameter Meanings
% data -> should be cells with a structure inside every cell that holds info for one Signal or Channel, might be single or multi-unit
% presentChoice -> Which choice of interest (1->R, 0->L)
% presentPDW -> High and Low plotting (0 -> No Just Choice, 1 -> Plot High and Low Bets)
% cohMarker - > Which coherence should we use 
% doWeNorm -> should we normalize data, mainly if we are pushing two sepearate Folder Files but that are of the same Day
% individualBlocks -> Does each data input belong to one individual Block,
% good for when blocks have two neurons in one channel
% plotFlag -> Should we Plot Spike Rate Comparison for Conditions
%% The data will be input as cells, and each cell will contain a structure
% (We have to somehow remove single PDW trials, might have to load PDS file
% for that)
%folderOfInterest = populationOne{7}; %Feed the Dates of interest
%%%%%%%%%%%% COMBINE FOLDERS IF SAME DAY LIKE YOU DID IN BEHAVIORRESULTS_PDW_CHOICES.MAT' %%%%%%%%%%%%%%%%%%%
%doWeNorm = 1;
allNeuronsSpikeCounts  = [];
allNeuronsDir = [];
channel = [];
selectivity={};
selectivity_ = [];
lastIndex = [];
index = 1;
for i = 1:length(data)
    motionStart = data{i}.openEvents.motionStart; %Select Period of Interest (Start)
    motionEnd = data{i}.openEvents.motionEnd; %Period of INterest (End)
    timeWindow = motionEnd-motionStart;
    spikeTimes = data{i}.spikeTimes; %This should be Nx1 (Spike Times)
    spikeCounts = sum(spikeTimes' > motionStart' & spikeTimes' <= motionEnd', 2)'; %Spike count in window of interest
    if wantFiringRates == 1
        spikeCounts = spikeCounts./timeWindow;
    end
    if doWeNorm == 1 % should we normalize?
        muSpikeCounts = mean(spikeCounts);
        nSpikeCounts = length(spikeCounts);
        steSpikeCounts = std(spikeCounts)/sqrt(nSpikeCounts);
        spikeCounts = (spikeCounts - muSpikeCounts)./steSpikeCounts;
    end
    if ~any(channel == data{i}.channel) || individualBlocks
        allNeuronsSpikeCounts(index,:) = spikeCounts; %Save all Neuron Spike Counts in Mat
        allNeuronsDir(index) = str2double(data{i}.dirSelectivity); %Their directions
        channel(index) = data{i}.channel;
        selectivity{index} = data{i}.dirSelectivity;
        selectivity_(index)  = str2double(data{i}.dirSelectivity);
        lastIndex(index) = length(spikeCounts);
        index = index + 1;
    else
        channelRow = channel == data{i}.channel;
        allNeuronsSpikeCounts(channelRow, lastIndex(channelRow)+1:lastIndex(channelRow)+length(spikeCounts)) = spikeCounts; %Save all Neuron Spike Counts in Mat
        lastIndex(channelRow) = lastIndex(channelRow)+length(spikeCounts);
    end
end
%% Eliminate any trials that one of the neurons does not exist, or better yet where one of the neurons has 0 counts
incompleteSpikeCount = any(allNeuronsSpikeCounts == 0, 1);
%% Plot Conditions for Spike Counts
% After Extracting Spike Counts Simply Align by COnditions of interest
% (which should be fed into function)
outlierMultiplier = 3; %Outlier multiplier (3 is standard)
% Simply Look at the behavior
choices = []; coherence=[]; direction=[]; signCoh=[]; pdw=[]; twoPDWs=[]; timeMotion=[];
index = 1;
allDates = {};
for n = 1:length(data) %Itereate through all Data we have
    if ~any(ismember(allDates, data{n}.date)) %If there is a match then do this 
        trialNums = length(data{n}.openEvents.choice);
        if max(data{n}.openEvents.choice) == 2
            choices(end+1:end+trialNums) = data{n}.openEvents.choice - 1;
        else
            choices(end+1:end+trialNums) = data{n}.openEvents.choice;
        end
        coherence(end+1:end+trialNums) = data{n}.openEvents.coherence;
        direction(end+1:end+trialNums) = data{n}.openEvents.direction;
        timeMotion(end+1:end+trialNums) = data{i}.openEvents.motionEnd - data{i}.openEvents.motionStart; %Time of Dot Presentation, or RT if RT setting
        signCoh(end+1:end+trialNums) = data{n}.openEvents.coherence .* cosd(data{n}.openEvents.direction);
        pdw(end+1:end+trialNums) = data{n}.openEvents.pdw;
        twoPDWs(end+1:end+trialNums) = data{n}.twoPDWs;
        allDates{index} = data{n}.date;
        index = index+1;
    end
end
%% Conditioning Spike Counts
%First Condition
%cohMarker = [0 .032 .064 .128 .256 .512]; %presentPDW = 0; %Present PDW data %presentChoice = 0; %Which choice to present
if presentPDW == 1 && presentChoice == 1
    choiceMarkerRight = 1;
    choiceMarkerLeft = 1;
    pdwMarkerHigh = 1;
    pdwMarkerLow = 0;
elseif presentPDW == 1 && presentChoice == 0
    choiceMarkerRight = 0;
    choiceMarkerLeft = 0;
    pdwMarkerHigh = 1;
    pdwMarkerLow = 0;
elseif presentPDW == 0 
    choiceMarkerRight = 1;
    choiceMarkerLeft = 0;
    pdwMarkerHigh = [0 1];
    pdwMarkerLow = [0 1];
end
condition = ismember(pdw, pdwMarkerHigh) & ismember(choices, choiceMarkerRight) & twoPDWs == 1 & ismember(coherence, cohMarker) & ~incompleteSpikeCount;
timeMotionRight = timeMotion(condition);
spikeCountRight = allNeuronsSpikeCounts(:, condition); %SpikeCounts for these conditions
%Find mean and Standard Deviation To remove outliers
muSpikeCountRight = mean(spikeCountRight,2); %mean
stdSpikeCountRight = std(spikeCountRight,[],2); %std
indexOutlier = any(spikeCountRight >= muSpikeCountRight + stdSpikeCountRight*outlierMultiplier |... 
    spikeCountRight <= muSpikeCountRight - stdSpikeCountRight*outlierMultiplier);%Any outliers? if so remove all Spikes for that trial
if any(indexOutlier)
    spikeCountRight(:, indexOutlier) = []; %New Samples
    muSpikeCountRight = mean(spikeCountRight, 2); %new Mean
    stdSpikeCountRight = std(spikeCountRight, [], 2); %new standard deviation
end
%Second Condition
condition = ismember(pdw, pdwMarkerLow) & ismember(choices, choiceMarkerLeft) & twoPDWs == 1 & ismember(coherence, cohMarker) & ~incompleteSpikeCount;
timeMotionLeft = timeMotion(condition);
spikeCountLeft = allNeuronsSpikeCounts(:, condition); %SpikeCounts for these contisions
%Find mean and Standard Deviation To remove outliers
muSpikeCountLeft = mean(spikeCountLeft,2); %mean
stdSpikeCountLeft = std(spikeCountLeft,[],2); %std
indexOutlier = any(spikeCountLeft >= muSpikeCountLeft + stdSpikeCountLeft*outlierMultiplier |... 
    spikeCountLeft <= muSpikeCountLeft - stdSpikeCountLeft*outlierMultiplier);%Any outliers? if so remove all Spikes for that trial
if any(indexOutlier)
    spikeCountLeft(:, indexOutlier) = []; %New Samples
    muSpikeCountLeft = mean(spikeCountLeft, 2); %new Mean
    stdSpikeCountLeft = std(spikeCountLeft, [], 2); %new standard deviation
end
%% Plot
S = 50;
%plotFlag = 1;
if size(spikeCountRight,1) == 3 && plotFlag == 1
    figure(1);
    scatter3(spikeCountRight(1,:), spikeCountRight(2,:), spikeCountRight(3,:), 'r'); hold on;
    scatter3(spikeCountLeft(1,:), spikeCountLeft(2,:), spikeCountLeft(3,:), 'b'); 
    scatter3(muSpikeCountRight(1), muSpikeCountRight(2), muSpikeCountRight(3), S, 'r', 'filled', 'MarkerEdgeColor', 'k'); 
    scatter3(muSpikeCountLeft(1), muSpikeCountLeft(2), muSpikeCountLeft(3), S, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold off;
    xlabel(selectivity{1}); ylabel(selectivity{2}); zlabel(selectivity{3});
elseif size(spikeCountRight,1) == 2 && plotFlag == 1
    figure(1);
    scatter(spikeCountRight(1,:), spikeCountRight(2,:), 'r'); hold on;
    scatter(spikeCountLeft(1,:), spikeCountLeft(2,:), 'b'); 
    scatter(muSpikeCountRight(1), muSpikeCountRight(2), S, 'r', 'filled', 'MarkerEdgeColor', 'k'); 
    scatter(muSpikeCountLeft(1), muSpikeCountLeft(2), S, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold off;
    xlabel(selectivity{1}); ylabel(selectivity{2}); 
elseif size(spikeCountRight,1) == 1 && plotFlag == 1
    figure(1);
    scatter(5, spikeCountRight(1,:), 'r'); hold on;
    scatter(5, spikeCountLeft(1,:), 'b'); 
    scatter(5, muSpikeCountRight(1), S, 'r', 'filled', 'MarkerEdgeColor', 'k'); 
    scatter(5, muSpikeCountLeft(1), S, 'b', 'filled', 'MarkerEdgeColor', 'k'); hold off;
    xlabel('Nothing'); ylabel(selectivity{1}); 
end
