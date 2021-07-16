%Miguel Vivar-Lazo
%3/2/2021
%Fetsch Lab: Making Pseudo Population

% 1) Keep Neurons together from Same Blocks
% 2) Normalize Firing Rates from all Neurons within Blocks
% 3) Partition Firing Rates by Coherence.*Direction and RT Quantiles
% 4) For each Neuron (within a block) combine the Firing Rates with other
% Neurons that fit the same conditions (Selectivity, Trials SignCoh and RT
% quantile and Choice/PDW)
% 5) Using these distributions for all possible combinations of
% dependencies fill out a PseudoPopulation Matrix (Neuron Selectivity -x-
% Trials). Try filling up as many trials with Firing Rates between Neurons
% that were recorded on the same session/block.
% 6) Fill in the holes in the pseudo populations by drawing from the
% distributions from all possible condition.
function [allNeuronsFiringRates, allNeuronsDir, direction, coherence, choices, reactionTime, pdw] = organizeForPseudoPop(data, doWeNorm, removeOutliers, individualBlocks)
%% 1) & 2)
%doWeNorm = 1;
allNeuronsFiringRates  = [];
allNeuronsDir = [];
channel = [];
selectivity={};
lastIndex = [];
index = 1;
for i = 1:length(data)
    motionStart = data{i}.openEvents.motionStart; %Select Period of Interest (Start)
    motionEnd = data{i}.openEvents.motionEnd; %Period of INterest (End)
    windowOfSpikes = motionEnd - motionStart;
    spikeTimes = data{i}.spikeTimes; %This should be Nx1 (Spike Times)
    spikeCounts = sum(spikeTimes' > motionStart' & spikeTimes' <= motionEnd', 2)'; %Spike count in window of interest
    spikeFR = spikeCounts./windowOfSpikes; 
    if doWeNorm == 1 % should we normalize?
        muSpikeFR = mean(spikeFR);
        nSpikeFR = length(spikeFR);
        steSpikeFR = std(spikeFR)/sqrt(nSpikeFR);
        spikeFR = (spikeFR - muSpikeFR)./steSpikeFR;
    end
    if ~any(channel == data{i}.channel) || individualBlocks
        allNeuronsFiringRates(index,:) = spikeFR; %Save all Neuron Spike Counts in Mat
        selectivity{index} = data{i}.dirSelectivity; %This line and the next one are the same 
        allNeuronsDir(index) = str2double(data{i}.dirSelectivity); %Their directions
        if isnan(allNeuronsDir(index))
            temp = data{i}.dirSelectivity;
            k = strfind(temp, '_');
            allNeuronsDir(index) = str2double(data{i}.dirSelectivity{1}(1:k{1}-1));
        end
        channel(index) = data{i}.channel;
        lastIndex(index) = length(spikeFR);
        index = index + 1;
    else
        channelRow = channel == data{i}.channel;
        allNeuronsFiringRates(channelRow, lastIndex(channelRow)+1:lastIndex(channelRow)+length(spikeCounts)) = spikeFR; %Save all Neuron Spike Counts in Mat
        lastIndex(channelRow) = lastIndex(channelRow)+length(spikeFR);
    end
end
%% Eliminate any trials that one of the neurons does not exist, or better yet where one of the neurons has 0 counts
incompleteSpikeCount = any(allNeuronsFiringRates == 0, 1);
%% Plot Conditions for Spike Counts
% After Extracting Spike Counts Simply Align by COnditions of interest
% (which should be fed into function)
% Simply Look at the behavior
choices = []; coherence=[]; direction=[]; signCoh=[]; pdw=[]; twoPDWs=[]; timeMotion=[];
index = 1;
allDates = {};
for n = 1:length(data) %Itereate through all Data we have
    if ~any(ismember(allDates, data{n}.date)) %If there is a match then do this 
        trialNums = length(data{n}.openEvents.choice);
        choices(end+1:end+trialNums) = data{n}.openEvents.choice - 1;
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
%% Remove any outliers and any Trials that dont have both PDWs
%removeOutliers = 1;
multiplierOutlier = 3;
if removeOutliers == 1
    muFiringRates = mean(allNeuronsFiringRates, 2);
    nFiringRates = size(allNeuronsFiringRates, 2);
    stdFiringRates = std(allNeuronsFiringRates, [], 2);
    outlierDetected = any(allNeuronsFiringRates > muFiringRates + stdFiringRates*multiplierOutlier | allNeuronsFiringRates < muFiringRates - stdFiringRates*multiplierOutlier, 1); %any trials with outliers in any neuron remove those
end
%% Return the results (Remove Outlier Trials & twoPDWs == 0 from Firing Rates, RT, Coherence, Choice, and PDW)
condition = ~outlierDetected & twoPDWs == 1;
allNeuronsFiringRates   = allNeuronsFiringRates(:,condition);
direction               = direction(condition);
coherence               = coherence(condition);
choices                 = choices(condition);
reactionTime            = timeMotion(condition);
pdw                     = pdw(condition);




