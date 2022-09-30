%Miguel Vivar-Lazo
%3/1/2021
%Fetsch Lab: Mapping Regression Weights for Various Dates 
%% Load Data

clear; close all
load DataSetForWeightsAnalysis

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
% % populationNonePDPLPS = {'hanzo20201103Dots1504.PDS', 'hanzo20201113Dots1422.PDS', 'hanzo20201113Dots1530.PDS', 'hanzo20201215Dots1413.PDS'};
% % populationOnePDPLPS = {'hanzo20201029Dots1523.PDS', 'hanzo20201117Dots1419.PDS', 'hanzo20201117Dots1535.PDS', 'hanzo20201123Dots1728.PDS', 'hanzo20201124Dots1614.PDS', 'hanzo20201201Dots1532.PDS', 'hanzo20201201Dots1644.PDS', 'hanzo20201208Dots1427.PDS'};
% % populationTwoPDPLPS = {'hanzo20201109Dots1521.PDS', 'hanzo20201110Dots1506.PDS', 'hanzo20201112Dots1643.PDS', 'hanzo20201119Dots1544.PDS', 'hanzo20201210Dots1633.PDS', 'hanzo20201211Dots1520.PDS', 'hanzo20201211Dots1413.PDS', 'hanzo20201217Dots1519.PDS', 'hanzo20201218Dots1348.PDS', 'hanzo20201218Dots1455.PDS'};
% % populationAllPLDSPS = [populationNonePDPLPS, populationOnePDPLPS, populationTwoPDPLPS];
% % %We will Look For Channels of that are spike sorted and contain 'exp' on
% % %them.
% % folderOfInterest = populationAll; %2 & 21 , no good
% % pldapsFolder = populationAllPLDSPS; %2
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
% %         if contains(filesInFolder(n).name, date) && contains(filesInFolder(n).name, '_exp') && contains(filesInFolder(n).name, '_Sorting_') && contains(filesInFolder(n).name, '_D')%Make sure File has the Same Data, 'Exp', 'Sorting', and Direction
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
% %                 %keyboard
% %                 if any(openEvents.blockNum > 1) 
% %                     if length(unique(openEvents.direction(openEvents.blockNum==2))) == 2 %No more than 2 directions
% %                         %keyboard %Make sure PLDAPS is same length as OEphys file
% %                         openEvents.motionStart      = openEvents.motionStart(openEvents.blockNum==2);
% %                         openEvents.motionEnd        = openEvents.motionEnd(openEvents.blockNum==2);
% %                         openEvents.FixPoint         = openEvents.FixPoint(openEvents.blockNum==2);
% %                         openEvents.FPAcquire        = openEvents.FPAcquire(openEvents.blockNum==2);
% %                         openEvents.endAcquire       = openEvents.endAcquire(openEvents.blockNum==2);
% %                         openEvents.speed            = openEvents.speed(openEvents.blockNum==2);
% %                         openEvents.diameter         = openEvents.diameter(openEvents.blockNum==2);
% %                         openEvents.trialNum         = openEvents.trialNum(openEvents.blockNum==2);
% %                         openEvents.direction        = openEvents.direction(openEvents.blockNum==2);
% %                         openEvents.timeDirection    = openEvents.timeDirection(openEvents.blockNum==2);
% %                         openEvents.trialStart       = openEvents.trialStart(openEvents.blockNum==2);
% %                         openEvents.coherence        = openEvents.coherence(openEvents.blockNum==2);
% %                         openEvents.choice           = openEvents.choice(openEvents.blockNum==2);
% %                         openEvents.pdw              = openEvents.pdw(openEvents.blockNum==2);
% %                         dataNeuron{index}.openEvents = openEvents; %replace it with new 
% %                     end
% %                 elseif length(tempFile.trialNum) > length(openEvents.motionStart)
% %                     tempFile.twotargconfidence = tempFile.twotargconfidence(ismember(tempFile.trialNum, openEvents.trialNum)); %Replace the twoConfidence vector by removing the values that belong to trials that openEvnts doesnt have                   
% %                 end      
% %             end
% %             dataNeuron{index}.twoPDWs = tempFile.twotargconfidence; %1=two targets, 0 = 1 target
% %             index = 1 + index;
% %         end
% %     end
% % end
%% Function for Spike Rates for Conditions Of Interest
presentChoice = 1; %Which Choice to concentrate on
presentPDW = 0; %Check PDW or check Choices
cohMarker = [0 .032 .064];
indAnalysis = 1;
doWeNorm = 1;
wantFiringRates = 1;
plotFlag = 0;
perQuant = .75;
weigths270Choice=[];weigths270PDWR=[];weigths270PDWL=[]; weights90Choice=[];weights90PDWR=[];weights90PDWL=[];
weights180Choice=[];weights180PDWR=[];weights180PDWL=[]; weights0Choice=[];weights0PDWR=[];weights0PDWL=[];
weigths225Choice=[];weigths225PDWR=[];weigths225PDWL=[]; weights135Choice=[];weights135PDWR=[];weights135PDWL=[];
weights45Choice=[];weights45PDWR=[];weights45PDWL=[]; weights315Choice=[];weights315PDWR=[];weights315PDWL=[];
rt270High=[];rt270Low=[];rt90High=[];rt90Low=[];rt0High=[];rt0Low=[];rt180High=[];rt180Low=[];
for i = [1, 3:10, 12:18, 20:length(folderOfInterest)] %7 has no Low Right CHoices
    dataFromSameDate = []; %Data cells with same Folder Date
    for n = 1:length(dataNeuron)
        if dataNeuron{n}.date == folderOfInterest{i}(1:end-4)
            dataFromSameDate(end+1) = n;
        end
    end
    
        data = dataNeuron(dataFromSameDate);
        [spikeCountHighRight, spikeCountLowRight, dirSelectivity, timeHighRight, timeLowRight] = MultiNeuron_TotalSpikeAnalysis(data, presentChoice, 1, cohMarker, doWeNorm, wantFiringRates, indAnalysis, plotFlag); %PDW + Right
        [spikeCountHighLeft, spikeCountLowLeft, dirSelectivity, timeHighLeft, timeLowLeft] = MultiNeuron_TotalSpikeAnalysis(data, 0, 1, cohMarker, doWeNorm, wantFiringRates, indAnalysis, plotFlag); %PDW + Left
        [spikeCountRight, spikeCountLeft, dirSelectivity, timeRight, timeLeft] = MultiNeuron_TotalSpikeAnalysis(data, presentChoice, presentPDW, cohMarker, doWeNorm, wantFiringRates, indAnalysis, plotFlag); %Choice
        
        % Logistic Regression (PDW Left)
        independentVar = [spikeCountRight spikeCountLeft]';
        dependentVar = [1*ones(1, size(spikeCountRight,2)) 2*ones(1, size(spikeCountLeft,2))]'; %I am pretty sure Right choices == 1 and Left == 2
        %[coefficientsChoice, dev, stats] = mnrfit(independentVar, dependentVar);
        [coefficientsChoice] = binarylogll_logit(independentVar, abs(dependentVar-2));
        % Logistic Regression (PDW Right)
        independentVar = [spikeCountHighRight spikeCountLowRight]';
        dependentVar = [1*ones(1, size(spikeCountHighRight,2)) 2*ones(1, size(spikeCountLowRight,2))]'; %I am pretty sure Right choices == 1 and Left == 2
        %[coefficientsPDWRight, dev, stats] = mnrfit(independentVar, dependentVar);
        [coefficientsPDWRight] = binarylogll_logit(independentVar, abs(dependentVar-2));
        % Logistic Regression (Choice)
        independentVar = [spikeCountHighLeft spikeCountLowLeft]';
        dependentVar = [1*ones(1, size(spikeCountHighLeft,2)) 2*ones(1, size(spikeCountLowLeft,2))]'; %I am pretty sure Right choices == 1 and Left == 2
        %[coefficientsPDWLeft, dev, stats] = mnrfit(independentVar, dependentVar);
        [coefficientsPDWLeft] = binarylogll_logit(independentVar, abs(dependentVar-2));

        %Grab the Weight for every Signal that is 270 or 90 and compare it to its weight from the Choice selection
        %try
        if any(dirSelectivity == 270)
            weigths270Choice(end+1:end+sum(dirSelectivity == 270)) = coefficientsChoice(find(dirSelectivity == 270)+1);
            weigths270PDWR(end+1:end+sum(dirSelectivity == 270)) = coefficientsPDWRight(find(dirSelectivity == 270)+1);
            weigths270PDWL(end+1:end+sum(dirSelectivity == 270)) = coefficientsPDWLeft(find(dirSelectivity == 270)+1);
            % Put Times (RT) in a vector %in order to map the time that
            % goes with each Weight
            spikeCount270 = [spikeCountRight(dirSelectivity == 270,:) spikeCountLeft(dirSelectivity == 270,:)]; %Spike Count for 270 Neuron
            rtChoice = [timeRight timeLeft];%Rt times align accordngly to what is above this line
            spCountIndecesForHighQuantile = find(any(spikeCount270 > quantile(spikeCount270, perQuant,2), 1)); %Index values of quantile above whatever percentage 
            spCountIndecesForLowQuantile = find(any(spikeCount270 > quantile(spikeCount270, 1-perQuant,2), 1)); %Index values of quantile above whatever percentage 
            %Find the highest quartile
            rt270High(end+1:end+length(rtChoice(spCountIndecesForHighQuantile))) = rtChoice(spCountIndecesForHighQuantile); %RT for these specific Quantiles of Spikes (high Spike Count)
            rt270Low(end+1:end+length(rtChoice(spCountIndecesForLowQuantile))) = rtChoice(spCountIndecesForLowQuantile); %RT for these specific Quantiles of Spikes (Low Spike Count)
        end
        if any(dirSelectivity == 90)
            weights90Choice(end+1:end+sum(dirSelectivity == 90)) = coefficientsChoice(find(dirSelectivity == 90)+1);
            weights90PDWR(end+1:end+sum(dirSelectivity == 90)) = coefficientsPDWRight(find(dirSelectivity == 90)+1);
            weights90PDWL(end+1:end+sum(dirSelectivity == 90)) = coefficientsPDWLeft(find(dirSelectivity == 90)+1);
            % Put Times (RT) in a vector
            spikeCount90 = [spikeCountRight(dirSelectivity == 90,:) spikeCountLeft(dirSelectivity == 90,:)]; %Spike Count for 270 Neuron
            rtChoice = [timeRight timeLeft];%Rt times align accordngly to what is above this line
            spCountIndecesForHighQuantile = find(any(spikeCount90 > quantile(spikeCount90, perQuant, 2), 1)); %Index values of quantile above whatever percentage 
            spCountIndecesForLowQuantile = find(any(spikeCount90 > quantile(spikeCount90, 1-perQuant, 2), 1)); %Index values of quantile above whatever percentage 
            %Find the highest quartile
            rt90High(end+1:end+length(rtChoice(spCountIndecesForHighQuantile))) = rtChoice(spCountIndecesForHighQuantile); %RT for these specific Quantiles of Spikes (high Spike Count)
            rt90Low(end+1:end+length(rtChoice(spCountIndecesForLowQuantile))) = rtChoice(spCountIndecesForLowQuantile); %RT for these specific Quantiles of Spikes (Low Spike Count)
        end
        if any(dirSelectivity == 0)
            weights0Choice(end+1:end+sum(dirSelectivity == 0)) = coefficientsChoice(find(dirSelectivity == 0)+1);
            weights0PDWR(end+1:end+sum(dirSelectivity == 0)) = coefficientsPDWRight(find(dirSelectivity == 0)+1);
            weights0PDWL(end+1:end+sum(dirSelectivity == 0)) = coefficientsPDWLeft(find(dirSelectivity == 0)+1);
            % Put Times (RT) in a vector
            spikeCount0 = [spikeCountRight(dirSelectivity == 0,:) spikeCountLeft(dirSelectivity == 0,:)]; %Spike Count for 270 Neuron
            rtChoice = [timeRight timeLeft];%Rt times align accordngly to what is above this line
            spCountIndecesForHighQuantile = find(any(spikeCount0 > quantile(spikeCount0, perQuant, 2), 1)); %Index values of quantile above whatever percentage 
            spCountIndecesForLowQuantile = find(any(spikeCount0 > quantile(spikeCount0, 1-perQuant, 2), 1)); %Index values of quantile above whatever percentage 
            %Find the highest quartile
            rt0High(end+1:end+length(rtChoice(spCountIndecesForHighQuantile))) = rtChoice(spCountIndecesForHighQuantile); %RT for these specific Quantiles of Spikes (high Spike Count)
            rt0Low(end+1:end+length(rtChoice(spCountIndecesForLowQuantile))) = rtChoice(spCountIndecesForLowQuantile); %RT for these specific Quantiles of Spikes (Low Spike Count)
        end
        if any(dirSelectivity == 180)
            weights180Choice(end+1:end+sum(dirSelectivity == 180)) = coefficientsChoice(find(dirSelectivity == 180)+1);
            weights180PDWR(end+1:end+sum(dirSelectivity == 180)) = coefficientsPDWRight(find(dirSelectivity == 180)+1);
            weights180PDWL(end+1:end+sum(dirSelectivity == 180)) = coefficientsPDWLeft(find(dirSelectivity == 180)+1);
            % Put Times (RT) in a vector
            spikeCount180 = [spikeCountRight(dirSelectivity == 180,:) spikeCountLeft(dirSelectivity == 180,:)]; %Spike Count for 270 Neuron
            rtChoice = [timeRight timeLeft];%Rt times align accordngly to what is above this line
            spCountIndecesForHighQuantile = find(any(spikeCount180 > quantile(spikeCount180, perQuant, 2),1)); %Index values of quantile above whatever percentage 
            spCountIndecesForLowQuantile = find(any(spikeCount180 > quantile(spikeCount180, 1-perQuant, 2), 1)); %Index values of quantile above whatever percentage 
            %Find the highest quartile
            rt180High(end+1:end+length(rtChoice(spCountIndecesForHighQuantile))) = rtChoice(spCountIndecesForHighQuantile); %RT for these specific Quantiles of Spikes (high Spike Count)
            rt180Low(end+1:end+length(rtChoice(spCountIndecesForLowQuantile))) = rtChoice(spCountIndecesForLowQuantile); %RT for these specific Quantiles of Spikes (Low Spike Count)
        end
        if any(dirSelectivity == 225)
            weigths225Choice(end+1:end+sum(dirSelectivity == 225)) = coefficientsChoice(find(dirSelectivity == 225)+1);
            weigths225PDWR(end+1:end+sum(dirSelectivity == 225)) = coefficientsPDWRight(find(dirSelectivity == 225)+1);
            weigths225PDWL(end+1:end+sum(dirSelectivity == 225)) = coefficientsPDWLeft(find(dirSelectivity == 225)+1);
        end
        if any(dirSelectivity == 135)
            weights135Choice(end+1:end+sum(dirSelectivity == 135)) = coefficientsChoice(find(dirSelectivity == 135)+1);
            weights135PDWR(end+1:end+sum(dirSelectivity == 135)) = coefficientsPDWRight(find(dirSelectivity == 135)+1);
            weights135PDWL(end+1:end+sum(dirSelectivity == 135)) = coefficientsPDWLeft(find(dirSelectivity == 135)+1);
        end
        if any(dirSelectivity == 45)
            weights45Choice(end+1:end+sum(dirSelectivity == 45)) = coefficientsChoice(find(dirSelectivity == 45)+1);
            weights45PDWR(end+1:end+sum(dirSelectivity == 45)) = coefficientsPDWRight(find(dirSelectivity == 45)+1);
            weights45PDWL(end+1:end+sum(dirSelectivity == 45)) = coefficientsPDWLeft(find(dirSelectivity == 45)+1);
        end
        if any(dirSelectivity == 315)
            weights315Choice(end+1:end+sum(dirSelectivity == 315)) = coefficientsChoice(find(dirSelectivity == 315)+1);
            weights315PDWR(end+1:end+sum(dirSelectivity == 315)) = coefficientsPDWRight(find(dirSelectivity == 315)+1);
            weights315PDWL(end+1:end+sum(dirSelectivity == 315)) = coefficientsPDWLeft(find(dirSelectivity == 315)+1);
        end
        %end
        % 
end
%[coefficients stats.p]
%% Plot (Do one for 0 and 180, angles that matter)
figure(10);
plot(([weigths270Choice]), ([weigths270PDWR]), 'bo'); hold on;
plot(([weigths270Choice]), ([weigths270PDWL]), 'b*');
plot(([weights90Choice]), ([weights90PDWR]), 'ro');
plot(([weights90Choice]), ([weights90PDWL]), 'r*');
plot(-.2:.01:.2, -.2:.01:.2, 'g--'); hold off
ylabel('Coefficient PDW')
xlabel('Coefficient Choice')
title('Comparison of Weights between Choice and PDW, In Orthogonal Units')
legend('270 Selective & Left PDW', '270 Selective & Right PDW', '90 Selective & Left PDW', '90 Selective & Right PDW')
grid off


figure(11);
plot(([weights0Choice]), ([weights0PDWL]), 'r*'); hold on;
plot(([weights0Choice]), ([weights0PDWR]), 'ro');
plot(([weights180Choice]), ([weights180PDWL]), 'bo');
plot(([weights180Choice]), ([ weights180PDWR]), 'b*');
plot(-.2:.01:.2, -.2:.01:.2, 'g--'); hold off
ylabel('Coefficient PDW')
xlabel('Coeffience Choice')
title('Comparison of Weights between Choice and PDW, In Parallel')
legend('0 Selective & Left PDW', '0 Selective & Preferred PDW', '180 Selective & Preferred PDW', '180 Selective & Right PDW')
grid off

% Add a linear regression line to this to see what the relationship between
% weights looks like?
weightsChoice = [weights0Choice weights180Choice weigths270Choice weights90Choice weights135Choice weigths225Choice weights45Choice weights315Choice];
weightsPDWR = [weights0PDWR weights180PDWR weigths270PDWR weights90PDWR weights135PDWR weigths225PDWR weights45PDWR weights315PDWR];
weightsPDWL = [weights0PDWL weights180PDWL weigths270PDWL weights90PDWL weights135PDWL weigths225PDWL weights45PDWL weights315PDWL];
betasLinear = fitlm(weightsChoice, weightsPDWR);
relationshipBtwChoicePDWRight = betasLinear.Coefficients.Estimate(2) .* (-0.2:.01:0.2);
betasLinear = fitlm(weightsChoice, weightsPDWL);
relationshipBtwChoicePDWLeft = betasLinear.Coefficients.Estimate(2) .* (-0.2:.01:0.2);
figure(12);
a = plot(weights0Choice, weights0PDWR , 'r.', 'Markersize' , 30); hold on;
b = plot(weights180Choice, weights180PDWR , 'b.', 'Markersize' , 30); 
c = plot(weigths270Choice, weigths270PDWR , 'g.', 'Markersize' , 30); 
m = plot(weights90Choice, weights90PDWR , 'm.', 'Markersize' , 30); 
k = plot(weights135Choice, weights135PDWR , 'k.', 'Markersize' , 30); 
y = plot(weigths225Choice, weigths225PDWR , 'y.', 'Markersize' , 30); 
c = plot(weights315Choice, weights315PDWR , 'c.', 'Markersize' , 30); 
w = plot(weights45Choice, weights45PDWR , '.', 'Color', [.5 .5 .5], 'Markersize' , 30); 
%Plot relationship
plot((-0.2:.01:0.2), relationshipBtwChoicePDWRight, 'r-');
plot(-.2:.01:.2, -.2:.01:.2, 'k--'); hold off
ylabel('Coefficient PDW')
xlabel('Coefficient Choice')
title('Choice vs Right PDW')
legend('0Deg', '180', '270', '90', '135', '225', '315', '45')
grid off

figure(13);
plot(weights0Choice, weights0PDWL , 'r*'); hold on;
plot(weights180Choice, weights180PDWL , 'b*');
plot(weigths270Choice, weigths270PDWL , 'g*'); 
plot(weights90Choice, weights90PDWL , 'm*');
plot(weights135Choice, weights135PDWL , 'k*'); 
plot(weigths225Choice, weigths225PDWL , 'y*');
plot(weights315Choice, weights315PDWL , 'c*'); 
plot(weights45Choice, weights45PDWL ,  '*', 'Color', [.5 .5 .5]); 
%Plot relationship
plot((-0.2:.01:0.2), relationshipBtwChoicePDWLeft, 'b-');
plot(-.2:.01:.2, -.2:.01:.2, 'k--'); hold off
title('Choice vs Left PDW')
% xlim([-0.2 0.2])
% ylim([-0.2 0.2])
ylabel('Coefficient PDW')
xlabel('Coefficient Choice')
legend('0Deg', '180', '270', '90', '135', '225', '315', '45')
grid off










