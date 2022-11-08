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
%% Load Data
% % nasOn = 1; %Are we looking through the NAS system
% % numPresentations = 1;
% % experiment = 1;
% % offlineSorting = 1;
% % spikeFile = 0;
% % contFile = 1;
% % %Folders of interest
% % populationNone = {'2020-11-03_15-04-31_exp', '2020-11-13_14-11-47_map', '2020-11-13_15-30-49_exp', '2020-12-15_14-13-19_exp'};
% % populationOne = {'2020-10-29_15-23-51_exp', '2020-11-17_14-19-33_exp', '2020-11-17_15-35-31_exp', '2020-11-23_17-29-03_exp', '2020-11-24_16-14-42_exp', '2020-12-01_15-32-32_exp', '2020-12-01_16-44-30_exp', '2020-12-08_14-27-29_exp'};
% % populationTwo = {'2020-11-09_15-22-01_exp', '2020-11-10_15-06-13_exp', '2020-11-12_16-43-23_exp', '2020-11-19_15-44-31_exp', '2020-12-10_16-28-50_exp', '2020-12-11_15-20-25_exp', '2020-12-11_13-45-56_map', '2020-12-17_15-19-51_exp', '2020-12-18_13-48-22_exp', '2020-12-18_14-56-00_exp'};
% % populationAll = [populationNone, populationOne, populationTwo];
% % %PLDAPS Files of Interest
% % % 'hanzo20201102Dots1409.PDS' + '2020-11-02_14-08-21_exp' file corrupted (pldaps one)
% % %'2020-11-12_14-48-49_exp' + 'hanzo20201112Dots_SALVAGED.mat  %fix the
% % %salvaged thing
% % populationNonePDPLPS = {'hanzo20201103Dots1504.PDS', 'hanzo20201113Dots1422.PDS', 'hanzo20201113Dots1530.PDS', 'hanzo20201215Dots1413.PDS'};
% % populationOnePDPLPS = {'hanzo20201029Dots1523.PDS', 'hanzo20201117Dots1419.PDS', 'hanzo20201117Dots1535.PDS', 'hanzo20201123Dots1728.PDS', 'hanzo20201124Dots1614.PDS', 'hanzo20201201Dots1532.PDS', 'hanzo20201201Dots1644.PDS', 'hanzo20201208Dots1427.PDS'};
% % populationTwoPDPLPS = {'hanzo20201109Dots1521.PDS', 'hanzo20201110Dots1506.PDS', 'hanzo202201112Dots1643.PDS', 'hanzo20201119Dots1544.PDS', 'hanzo20201210Dots1633.PDS', 'hanzo20201211Dots1520.PDS', 'hanzo20201211Dots1413.PDS', 'hanzo20201217Dots1519.PDS', 'hanzo20201218Dots1348.PDS', 'hanzo20201218Dots1455.PDS'};
% % populationAllPLDSPS = [populationNonePDPLPS, populationOnePDPLPS, populationTwoPDPLPS];
% % %We will Look For Channels of that are spike sorted and contain 'exp' on
% % %them.
% % folderOfInterest = populationOne;
% % pldapsFolder = populationOnePDPLPS;
% % % Search in FileName for 
% % dataNeuron = cell(1,1);
% % index = 1;
% % for i = 1:length(folderOfInterest)
% %     %Which directory
% %     if nasOn == 0
% %         dirname = ['C:\Users\fetschlab\data\' folderOfInterest{i} '\'];
% %         pldapsdir = 'Z:\data\hanzo\';
% %     else
% %         dirname = ['Z:\data\hanzo_neuro\' folderOfInterest{i} '\'];
% %         pldapsdir = 'Z:\data\hanzo\';
% %     end
% %     filesInFolder = dir(dirname); %Go to directory of Interest
% %     date = folderOfInterest{i}(1:end-4); %Remove last 4 char from str of folder (should either be '_map' or '_exp'
% %     for n = 1:length(filesInFolder)
% %         if contains(filesInFolder(n).name, date) && contains(filesInFolder(n).name, '_exp') && contains(filesInFolder(n).name, '_Sorting_')%Make sure File has the Same Data, 'Exp', 'Sorting', and Direction
% %             dataNeuron{index}.date = date; % date of the file
% %             dataNeuron{index}.nameFile = filesInFolder(n).name;%Take the Folder Name and Its Direction
% %             dataNeuron{index}.dirSelectivity = extractBetween(filesInFolder(n).name, "_D", ".psort"); %Pickout Direction
% %             dataNeuron{index}.channel = str2double(extractBetween(filesInFolder(n).name, "_Ch", "_Sorting"));
% %             %Load information
% %             pSorterResults = Psort_read_psort([dirname filesInFolder(n).name]); % Load Spike Files
% %             if isnan(str2double(extractBetween(filesInFolder(n).name, '_exp', '_Ch'))) %Which events file too load
% %                 jobID = 1;
% %             else 
% %                 jobID = str2double(extractBetween(filesInFolder(n).name, '_exp', '_Ch'));
% %             end
% %             [openEvents, ~, openCont, ~] = neuralExtractorMapping(jobID, numPresentations, dataNeuron{index}.channel, experiment, offlineSorting, spikeFile, contFile, [], dirname);
% %             %Save information
% %             dataNeuron{index}.openEvents = openEvents; %transfer name
% %             dataNeuron{index}.spikeTimes = openCont(pSorterResults.topLevel_data.ss_index == 1); %spike Times from offline sorter   
% %             %Load PDS file
% %             tempFile = SaccadeTraining_BehaviorResults(pldapsdir, pldapsFolder{i}(1:5), str2double(pldapsFolder{i}(6:13)), pldapsFolder{i}(14:17), pldapsFolder{i}(18:21)); %You need to automate this
% %             if length(tempFile.twotargconfidence) ~= length(openEvents.motionStart)
% %                 if any(openEvents.blockNum > 1) 
% %                     %keyboard %Make sure PLDAPS is same length as OEphys file
% %                     openEvents.motionStart      = openEvents.motionStart(openEvents.blockNum==2);
% %                     openEvents.motionEnd        = openEvents.motionEnd(openEvents.blockNum==2);
% %                     openEvents.FixPoint         = openEvents.FixPoint(openEvents.blockNum==2);
% %                     openEvents.FPAcquire        = openEvents.FPAcquire(openEvents.blockNum==2);
% %                     openEvents.endAcquire       = openEvents.endAcquire(openEvents.blockNum==2);
% %                     openEvents.speed            = openEvents.speed(openEvents.blockNum==2);
% %                     openEvents.diameter         = openEvents.diameter(openEvents.blockNum==2);
% %                     openEvents.trialNum         = openEvents.trialNum(openEvents.blockNum==2);
% %                     openEvents.direction        = openEvents.direction(openEvents.blockNum==2);
% %                     openEvents.timeDirection    = openEvents.timeDirection(openEvents.blockNum==2);
% %                     openEvents.trialStart       = openEvents.trialStart(openEvents.blockNum==2);
% %                     openEvents.coherence        = openEvents.coherence(openEvents.blockNum==2);
% %                     openEvents.choice           = openEvents.choice(openEvents.blockNum==2);
% %                     openEvents.pdw              = openEvents.pdw(openEvents.blockNum==2);
% %                     dataNeuron{index}.openEvents = openEvents; %replace it with new 
% %                 elseif length(tempFile.trialNum) > length(openEvents.motionStart)
% %                     tempFile.twotargconfidence = tempFile.twotargconfidence(ismember(tempFile.trialNum, openEvents.trialNum)); %Replace the twoConfidence vector by removing the values that belong to trials that openEvnts doesnt have                   
% %                 end               
% %             end
% %             dataNeuron{index}.twoPDWs = tempFile.twotargconfidence; %1=two targets, 0 = 1 target
% %             index = 1 + index;
% %         end
% %     end
% % end
%% Make Variables for all possible distributions of FR (This if from PLDAPS algorithm in setup, good thing to know)
cc.rt = [1, 2, 3];
cc.Coherences = [0, 0.032, 0.064 0.128 .256 .512]; % 0.128, 0.256, 0.512];
cc.Direction = [0, 180];
cc.Choice = [0, 1];
cc.PDW = [0, 1];
cc.Selectivity = [0, 45, 90, 135, 180, 270, 315]; %225

fn = fieldnames(cc);
numConds = 1;
for parameter=1:length(fn)
    numConds = numConds*length(cc.(fn{parameter}));
end
c=cell(1,numConds);
numCondsTillNow=1;
for parameter=1:length(fn)
    numParmValues=length(cc.(fn{parameter}));
    for condition = 1:numConds
        thisParmValue=floor(mod((condition-1)/numCondsTillNow,numParmValues)+1);
        if(isnumeric(cc.(fn{parameter})) || islogical(cc.(fn{parameter})))
            c{condition}.(fn{parameter})=cc.(fn{parameter})(thisParmValue);
        else
            fn2=fieldnames(cc.(fn{parameter}){thisParmValue});
            for ParmValueField=1:length(fn2)
                c{condition}.(fn2{ParmValueField})=cc.(fn{parameter}){thisParmValue}.(fn2{ParmValueField});
            end
        end
        c{condition}.FiringRate = []; %Blanks for distribution in a population
    end
    numCondsTillNow = numCondsTillNow*numParmValues;
end
%% Organize it (Gotta run this one and the block before or you will get an error)
doWeNorm = 1;
removeOutliers = 1;
rtQuantiles = length(cc.rt) - 1;
individualBlocks = 1;
numTotalTrials = 12000;
PopulationMatrix = nan(length(cc.Selectivity), numTotalTrials); %keeping track of the conditions for each pseudo trial
PopulationMatrixRT = nan(1, numTotalTrials);
PopulationMatrixCoherence = nan(1, numTotalTrials);
PopulationMatrixDirection = nan(1, numTotalTrials);
PopulationMatrixChoice = nan(1, numTotalTrials);
PopulationMatrixPDW = nan(1, numTotalTrials);
index = 1; %keeping track of the latest index for population matrix
for i = [1, 3:10, 12:18, 20:length(folderOfInterest)] %1:length(folderOfInterest)
    dataFromSameDate = []; %Data cells with same Folder Date
    for n = 1:length(dataNeuron) %Find all files with the same data, that way you can combine them because they will probably be recorded together
        if dataNeuron{n}.date == folderOfInterest{i}(1:end-4)
            dataFromSameDate(end+1) = n;
        end
    end
    % Organize data, using the function we used
    [allNeuronsFiringRates, allNeuronsDir, direction, coherence, choice, reactionTime, pdw] = organizeForPseudoPop(dataNeuron(dataFromSameDate), doWeNorm, removeOutliers, individualBlocks);
    
    %Break RTs, into quantiles
    rtPartitionMarkers = quantile(reactionTime, rtQuantiles); 
    rtGroup = zeros(size(reactionTime)); %Make a vector for RT categories
    for ii = 1:length(rtPartitionMarkers)
        if ii == 1
            rtGroup(reactionTime <= rtPartitionMarkers(ii)) = ii;
        elseif ii == length(rtPartitionMarkers) %last one
            rtGroup(reactionTime > rtPartitionMarkers(ii)) = ii+1;
            rtGroup(reactionTime > rtPartitionMarkers(ii-1) & reactionTime <= rtPartitionMarkers(ii)) = ii;
        else %anything in between
            rtGroup(reactionTime > rtPartitionMarkers(ii-1) & reactionTime <= rtPartitionMarkers(ii)) = ii;
        end
    end
    %Group Data
    tic
    uniqDirection = unique(direction); uniqCoh = unique(coherence); uniqChoice = unique(choice); uniqRT = unique(rtGroup); uniqPDW = unique(pdw);
    nTrialsForBlock = length(direction); %number of trials in block
    for d = 1:length(allNeuronsDir)
        selectivityOfInterest = allNeuronsDir(d);
        if 0 || selectivityOfInterest == 225
            break 
        end
        for z = uniqDirection %Iterate through all trials, direction should be the length of trials and it should be the same length as the rest of PDW, coherence, and so on
            for y = uniqCoh
                for x = uniqChoice
                    for w = uniqRT
                        for u = uniqPDW
                            for n = 1:length(c) %Finding Group with distribution of interest
                                if c{n}.rt == w && c{n}.Coherences == y && c{n}.Direction == z && c{n}.Choice == x && c{n}.PDW == u && c{n}.Selectivity == selectivityOfInterest
                                    condition = direction == z & coherence == y & choice == x & rtGroup == w & pdw == u; %Find trials with this condition
                                    c{n}.FiringRate(end+1: end+sum(condition)) = allNeuronsFiringRates(d, condition); %Place firing rate in respective distribution
                                end
                            end
                        end
                    end
                end
            end
        end
        %After making seperate populations, make sure to now Save Firing Rates
        %in Blocks that have more than 2 neurons. In the big Matrix, there will
        %be some filled in and a lot of empties, make sure to fill in the
        %empties
        if length(allNeuronsDir) > 1 %If this block is recording from more than 1 neuron them lets put it in the population
            indexOfSelectivity = cc.Selectivity == selectivityOfInterest;
            PopulationMatrix(indexOfSelectivity, index : (index-1)+nTrialsForBlock) = allNeuronsFiringRates(d, :);
        end
    end
    toc
    %Also save one instance of the trial information for this 
    if length(allNeuronsDir) > 1 %If this block is recording from more than 1 neuron them lets put it in the population
        PopulationMatrixRT(index : (index-1)+nTrialsForBlock) = rtGroup; %Fix the way we match up the index
        PopulationMatrixCoherence(index : (index-1)+nTrialsForBlock) = coherence;
        PopulationMatrixDirection(index : (index-1)+nTrialsForBlock) = direction;
        PopulationMatrixChoice(index : (index-1)+nTrialsForBlock) = choice;
        PopulationMatrixPDW(index : (index-1)+nTrialsForBlock) = pdw;
        index = index + nTrialsForBlock + 1; %Next blank spot
    end

end
%% Remove any NaNs from vector (list of Coh, choice and so on)
removeNaNs = isnan(PopulationMatrixCoherence);
PopulationMatrixRT          = PopulationMatrixRT(~removeNaNs);
PopulationMatrixDirection   = PopulationMatrixDirection(~removeNaNs);
PopulationMatrixCoherence   = PopulationMatrixCoherence(~removeNaNs);
PopulationMatrixChoice      = PopulationMatrixChoice(~removeNaNs);
PopulationMatrixPDW         = PopulationMatrixPDW(~removeNaNs);
PopulationMatrix            = PopulationMatrix(:, ~removeNaNs);
%% Make new matrix variables that will be constantly replace when looping
tempPopulationMatrixRT          = PopulationMatrixRT;
tempPopulationMatrixDirection   = PopulationMatrixDirection;
tempPopulationMatrixCoherence   = PopulationMatrixCoherence;
tempPopulationMatrixChoice      = PopulationMatrixChoice;
tempPopulationMatrixPDW         = PopulationMatrixPDW;
tempPopulationMatrix            = PopulationMatrix;
%% After Making Distributions make Sure to fill in what is needed from the Blanks (This is what you want to run in a loop)
numBoot = 100; %normally 1 if you want to run it only once
plotflag = 0;
bootWeights_ = zeros(numBoot, length(cc.Selectivity)); bootWeightsHigh_ = zeros(numBoot, length(cc.Selectivity)); bootWeightsLow_ = zeros(numBoot, length(cc.Selectivity));
lookingForCoh = [0 .032 .064 .128 .256 .512]; %Which coherences are you looking for
for nb = 1:numBoot
    for d = 1:length(PopulationMatrixCoherence) %Iterate through all trials, direction should be the length of trials and it should be the same length as the rest of PDW, coherence, and so on
        interestedRT        = PopulationMatrixRT(d);
        interestedDir       = PopulationMatrixDirection(d);
        interestedCoh       = PopulationMatrixCoherence(d);
        interestedChoice    = PopulationMatrixChoice(d);
        interestedPDW       = PopulationMatrixPDW(d);
        emptyFiringRates = find(isnan(PopulationMatrix(:,d)))';
        for z = emptyFiringRates %Then iterate through all missing Firing Rate inputs
            whichSelectivity = cc.Selectivity(z); %Which selectivity direction is missing
            for n = 1:length(c) %Then iterate to find and input a psuedo Firing Rate from distributions
                if c{n}.rt == interestedRT && c{n}.Coherences == interestedCoh && c{n}.Direction == interestedDir && c{n}.Choice == interestedChoice && c{n}.PDW == interestedPDW && c{n}.Selectivity == whichSelectivity
                    try %if this doesnt work then its because there was no FR for that particular parameters for the other neuron, which is sad (but okay)
                        tempPopulationMatrix(z, d) =  c{n}.FiringRate(randi(length(c{n}.FiringRate))); %replace the blank FR for that Trial and that Selectivity with one randomly from the distribution I constructed
                    catch
                        %keyboard
                    end
                end
            end
        end
    end
% Remove any blanks, becasue their might be some
    areThereNaNsInPopulation = isnan(mean(tempPopulationMatrix)); %There will be some NaNs in smaller data sets
    %remove these
    tempPopulationMatrixRT          = PopulationMatrixRT(~areThereNaNsInPopulation);
    tempPopulationMatrixDirection   = PopulationMatrixDirection(~areThereNaNsInPopulation);
    tempPopulationMatrixCoherence   = PopulationMatrixCoherence(~areThereNaNsInPopulation);
    tempPopulationMatrixChoice      = PopulationMatrixChoice(~areThereNaNsInPopulation);
    tempPopulationMatrixPDW         = PopulationMatrixPDW(~areThereNaNsInPopulation);
    tempPopulationMatrix            = tempPopulationMatrix(:, ~areThereNaNsInPopulation);
% Now that you have this Population Matrix you can apply whatever you want to apply 
    condition_forChoice = ismember(tempPopulationMatrixCoherence, lookingForCoh);
    condition_forChoice_High = ismember(tempPopulationMatrixCoherence, lookingForCoh) & ismember(tempPopulationMatrixPDW, 1);
    condition_forChoice_Low = ismember(tempPopulationMatrixCoherence, lookingForCoh) & ismember(tempPopulationMatrixPDW, 0);
    %[betas_test, ~, stats] = mnrfit(PopulationMatrix(:,condition_forChoice)', abs(PopulationMatrixChoice(condition_forChoice)'+1)); % you want choices for right == 1 and left == 2%
    [betas_test] = binarylogll_logit(tempPopulationMatrix(:,condition_forChoice)', tempPopulationMatrixChoice(condition_forChoice)'); % you want choices for right == 1 and left == 2
    [betas_test_High] = binarylogll_logit(tempPopulationMatrix(:,condition_forChoice_High)', tempPopulationMatrixChoice(condition_forChoice_High)'); %R vs L on High Bets
    [betas_test_Low] = binarylogll_logit(tempPopulationMatrix(:,condition_forChoice_Low)', tempPopulationMatrixChoice(condition_forChoice_Low)'); %R vs L on Low Bets
    betas_test = betas_test(2:end)'; %remove the bias (for now)
    betas_test_High = betas_test_High(2:end)'; %remove the bias (for now)
    betas_test_Low = betas_test_Low(2:end)'; %remove the bias (for now)
    if plotflag == 1
        figure(91);
        plot([-fliplr(cc.Selectivity) cc.Selectivity(2:end)], [flip(betas_test) betas_test(2:end)], 'ro--');
    end
    % Place numerous Beta (weights) value iterations in a big matrix
    % (Iteration -x- Weights)
    bootWeights(nb, :) = betas_test; 
    bootWeightsHigh(nb, :) = betas_test_High;
    bootWeightsLow(nb, :) = betas_test_Low;
    fprintf('Iteration: %2.f\n', nb);
end
%% Plot the 'bootstrap' weights of the population
muBootWeights = mean(bootWeights);
stdBootWeights = std(bootWeights);
numBootWeights = size(bootWeights, 1);
steBootWeights = stdBootWeights./sqrt(numBootWeights);
figure(91);
errorbar([-fliplr(cc.Selectivity) cc.Selectivity(2:end)], [flip(muBootWeights) muBootWeights(2:end)], [flip(steBootWeights) steBootWeights(2:end)], 'bo--'); %hold off
title('Weight Coefficients for Population: On Choices')
xlabel('Motion Direction Selectivity')
ylabel('Weight Values (a.u)')
legend('All Coherence', 'Low Coherence')
%% Plot the 'bootstrap' weights of the population (High vs Low Bets Situation)
muBootWeightsHigh = mean(bootWeightsHigh);
numBootWeightsHigh = size(bootWeightsHigh, 1);
steBootWeightsHigh = muBootWeightsHigh./sqrt(numBootWeightsHigh);
muBootWeightsLow = mean(bootWeightsLow);
numBootWeightsLow = size(bootWeightsLow, 1);
steBootWeightsLow = muBootWeightsLow./sqrt(numBootWeightsLow);
figure(92);
errorbar([-fliplr(cc.Selectivity) cc.Selectivity(2:end)], [flip(muBootWeightsHigh) muBootWeightsHigh(2:end)], [flip(steBootWeightsHigh) steBootWeightsHigh(2:end)], 'rs-'); hold on 
errorbar([-fliplr(cc.Selectivity) cc.Selectivity(2:end)], [flip(muBootWeightsLow) muBootWeightsLow(2:end)], [flip(steBootWeightsLow) steBootWeightsLow(2:end)], 'bs-'); hold off
title('Weight Coefficients for Population: On Choices')
xlabel('Motion Direction Selectivity')
ylabel('Weight Values (a.u)')
legend('High Bet (Low Coh)', 'Low Bet (Low Coh)', 'High Bet (All Coh)', 'Low Bet (All Coh)')
%% Now do one for High and Low Confidence
% % % condition_forBetas = PopulationMatrixPDW == 1 & ismember(PopulationMatrixCoherence, [0 .032 .064]);
% % % betas_test_high = mnrfit(PopulationMatrix(:, condition_forBetas)', abs(PopulationMatrixPDW(condition_forBetas)'+1)); % you want choices for right == 1 and left == 2
% % % condition_forBetas = PopulationMatrixPDW == 0 & ismember(PopulationMatrixCoherence, [0 .032 .064]);
% % % betas_test_low = mnrfit(PopulationMatrix(:, condition_forBetas)', abs(PopulationMatrixChoice(condition_forBetas)'+1)); % you want choices for right == 1 and left == 2
% % % betas_test_high = betas_test_high(2:end)'; %remove the bias (for now)
% % % betas_test_low = betas_test_low(2:end)'; %remove the bias (for now)
% % % figure(92);
% % % plot([-fliplr(cc.Selectivity) cc.Selectivity(2:end)], [flip(betas_test_high) betas_test_high(2:end)], 'ro--'); hold on;
% % % plot([-fliplr(cc.Selectivity) cc.Selectivity(2:end)], [flip(betas_test_low) betas_test_low(2:end)], 'bo--'); hold off;
% % % legend('High', 'Low')





