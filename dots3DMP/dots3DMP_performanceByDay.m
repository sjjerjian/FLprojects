% 

% Not all of the days in this dataset have the cue-conflict condition, but
% keeping them in because we can still look at choice biases and thresholds
% over days

clear; clc

cd /Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data
load lucio_20210315-20210801_clean.mat
alldata = data;

conftask = 2; % 1=colorbars, 2=PDW
RTtask = 1;

%%

mods   = unique(data.modality); 
cohs   = unique(data.coherence); 
deltas = unique(data.delta);
hdgs   = unique(data.heading);

days = unique(data.date); %date
nDays = length(days);

%%

fnames = fieldnames(alldata);

% loop over days

clear nCorrect* nTrials*
clear parsedData gfit wves deltas_day
for d = 1:nDays
    
    clear data;
    trials_day = alldata.date == days(d);
    
    % cutoff for the min number of trials to get reliable estimates
    
%     if sum(trials_day) > 600
        
        for F = 1:length(fnames)
            data.(fnames{F}) = alldata.(fnames{F})(trials_day);
        end
        nTrials(d,1) = sum(trials_day);
        nCorrect(d,1) = sum(data.correct);
        
        
        nTrialsConf(d,1) = sum(data.PDW==1);
        nTrialsConf(d,2) = sum(data.PDW==0);
        
        nCorrectConf(d,1) = sum(data.correct & data.PDW==1);
        nCorrectConf(d,2) = sum(data.correct & data.PDW==0);
        
        deltas_day{d,1} = unique(data.delta);
        hasDelta(d) = any(data.delta~=0);
        
        parsedData(d) = dots3DMP_parseData(data,mods,cohs,deltas_day{d},hdgs,conftask,RTtask); % basic parsing of data, logistic fits
        
        gfit(d) = dots3DMP_fit_cgauss(data,mods,cohs,deltas_day{d},conftask,RTtask); % gaussian fitting
        
        wves(d) = dots3DMP_cueWeights(gfit(d),cohs,deltas_day{d});
       
            
%     end

end

%% plot example day


day = 38; % 
dots3DMP_plots_cgauss_func(gfit(day),parsedData(day),mods,cohs,deltas,hdgs,conftask,RTtask)

%% summary statistics over days

% loop over days to extract necessary values into a simpler format for
% plotting

accuracy = nCorrect./nTrials;
accuracyConf = nCorrectConf./nTrialsConf;

[wvesPred,wvesEmp] = deal(nan(nDays,length(cohs)));

[choiceBias,confBias,RTBias] = deal(nan(length(mods),length(cohs),length(deltas)+1,nDays));

clear *Bias

for d=1:nDays
    
    % thresholds

    
    
    % biases/shifts
    if length(deltas_day{d})>1
        choiceBias(:,:,1:length(deltas)+1,d) = gfit(d).choice.mu;
        confBias(:,:,1:length(deltas)+1,d)   = gfit(d).conf.mu;
        RTBias(:,:,1:length(deltas)+1,d)     = gfit(d).RT.mu;
    else
        choiceBias(:,:,length(deltas)+1,d)   = gfit(d).choice.mu(:,:,end);
        confBias(:,:,length(deltas)+1,d)   = gfit(d).conf.mu(:,:,end);
        RTBias(:,:,length(deltas)+1,d)   = gfit(d).RT.mu(:,:,end);
    end
    
    % wVes
    wvesPred(d,:) = wves(d).choice.pred;

    if length(deltas_day{d})>1
        wvesEmp(d,:)  = wves(d).choice.emp;

    end    
end


%%

% ACCURACY
figure; subplot(211); hold on
plot(1:nDays,accuracy,'k',1:nDays,accuracyConf(:,1),'r',1:nDays,accuracyConf(:,2),'r--')
legend('Mean acc','High Bet Acc','Low Bet Acc')
xlabel('Session no.'); ylabel('Accuracy');
subplot(212); hold on
hr=refline(1,0); hr.Color = 'r';
scatter(accuracyConf(:,2),accuracyConf(:,1));
axis([0.5 0.9 0.5 0.9])
xlabel('Low Bet Acc'); ylabel('High Bet Acc');
%%
% VESTIBULAR WEIGHT
figure; hold on;
plot(1:nDays,wvesPred(:,1),'k--',1:nDays,wvesPred(:,2),'k-')
plot(1:nDays,wvesEmp(:,1),'r--',1:nDays,wvesEmp(:,2),'r-')


%% BIASES

c = 2; % coh

figure; 
for m=1:length(mods)
    subplot(3,1,m); hold on
    semilogy(1:nDays,squeeze(choiceBias(m,c,4,:)),'k','linew',1.5);
    semilogy(1:nDays,squeeze(confBias(m,c,4,:)),'r','linew',1.5);
    semilogy(1:nDays,squeeze(RTBias(m,c,4,:)),'b','linew',1.5);
%     ylim([-5 5]);
end
legend('Choice','Conf','RT')

if 0
corrcoef(squeeze(choiceBias(m,c,4,:)),squeeze(confBias(m,c,4,:)))
corrcoef(squeeze(choiceBias(m,c,4,:)),squeeze(RTBias(m,c,4,:)))
corrcoef(squeeze(RTBias(m,c,4,:)),squeeze(confBias(m,c,4,:)))
end

figure; hold on;
scatter(squeeze(choiceBias(m,c,4,:)),squeeze(confBias(m,c,4,:)),100,'k','.')
lsline;
% scatter(squeeze(choiceBias(m,c,4,:)),squeeze(RTBias(m,c,4,:)),30,'r','.')
% scatter(squeeze(RTBias(m,c,4,:)),squeeze(confBias(m,c,4,:)),30,'b','.')

%% SHIFTS UNDER CONFLICT CONDITION



