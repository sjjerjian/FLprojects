% first-pass at reading in and plotting data from Marton et al.,
% "Validating a dimension of doubt in decision-making: A proposed
% endophenotype for obsessive-compulsive disorder"
% CF & KT 09-2021

clear
close all

% change this to the folder where you have the data/code
% cd '/MATLAB Drive'
cd '/Users/chris/Documents/MATLAB/Projects/DoubtConf(Nestadt)'
data = readtable('RDMPilotRDM_GIANT_VC_Repository.xlsx');
data = table2struct(data,'ToScalar',true);

% use 'unique' on the independent variables, both to make sure we
% understand what's in the data, and to use as loop index variables below
uID = unique(data.ID); % Subject Identification #
uAGELT50 = unique(data.AGELT50); %Age less than 50 years: 1= age>50, 2=age<50
uGroupx2 = unique(data.Groupx2); % Disorder group: 1 = control, 2 = OCD 
uTest = unique(data.Test); % 0=practice, 1=$$feedback, 2=confidence, 3=$$penalty
uCoh = unique(data.Coherence); % motion strength (% coherence)

goodTrials = ~isnan(data.Correct);
badTrials_conf = isnan(data.Confidence) & strcmp(data.Test,'test2');
goodTrials(badTrials_conf) = false;


%% summarize choice, RT, conf as a function of coherence, for the whole dataset

% preallocate
pctCor = nan(length(uCoh),1);
RTmean = nan(length(uCoh),1);
RTse = nan(length(uCoh),1);
confMean = nan(length(uCoh),1);
confSE = nan(length(uCoh),1);
N = nan(length(uCoh),1);

for c = 1:length(uCoh)
    
    % start w basic plots of %corr, RT, and conf as a function of Coherence
    I = goodTrials & data.Coherence==uCoh(c);
    N(c) = sum(I);
    pctCor(c) = sum(data.Correct(I)==2) / N(c); % 1=err, 2=corr
    
    RTmean(c) = mean(data.ReactionTime(I));
    RTse(c) = std(data.ReactionTime(I))/sqrt(N(c)); % standard approximation of standard error of the mean (SEM),
                                                    % used for plotting error bars and some statistical tests
    J = I & strcmp(data.Test,'test2');
    confMean(c) = nanmean(data.Confidence(J));
    confSE(c) = nanstd(data.Confidence(J))/sqrt(sum(J));

end
% standard error for a percentage has a convenient formula
% and can be done outside the loop:
% plotting pctCorSE with error bar
pctCorSE = sqrt( (pctCor.*(1-pctCor)) ./ N );

% plot it!
figure(1); set(gcf,'Color',[1 1 1],'Position',[400 700 900 300],'PaperPositionMode','auto'); clf;
subplot(1,3,1);
errorbar(uCoh,RTmean,RTse,'bo-')
    ylabel('Mean RT (ms)')
    xlabel('Motion Coherence (%)')
    title ('RT')

subplot(1,3,2);
errorbar(uCoh,pctCor,pctCorSE,'ro-')
    ylabel('Percent correct')
    xlabel('Motion Coherence (%)')
    title ('accuracy')

subplot(1,3,3);
K = ~isnan(confMean); % some cohs were not used in conf experiment and thus gave NaN in above loop; omit those
errorbar(uCoh(K),confMean(K),confSE(K),'go-')
    ylabel('Mean confidence report')
    xlabel('Motion Coherence (%)')
    title ('confidence')
     
    
    
%% repeat this but split by participant

% preallocate: now these will be 2D arrays (matrices): uCoh and uID
pctCor = nan(length(uCoh),length(uID));
RTmean = nan(length(uCoh),length(uID));
RTse = nan(length(uCoh),length(uID));
confMean = nan(length(uCoh),length(uID));
confSE = nan(length(uCoh),length(uID));
N = nan(length(uCoh),length(uID));


% for n = 1:length(uID)
for n = 1:8 % just look at a few participants for now, 81 is a lot! 
    
    for c = 1:length(uCoh)
        I = goodTrials & data.Coherence==uCoh(c) & data.ID==uID(n);
        N(c,n) = sum(I);
        
        pctCor(c,n) = sum(data.Correct(I)==2) / N(c,n); % 1=err, 2=corr

        RTmean(c,n) = mean(data.ReactionTime(I));
        RTse(c,n) = std(data.ReactionTime(I))/sqrt(N(c,n));
                                                          
        J = I & strcmp(data.Test,'test2');
        if sum(J)>1
            confMean(c,n) = nanmean(data.Confidence(J));
            confSE(c,n) = nanstd(data.Confidence(J))/sqrt(sum(J));
        end
    end

    pctCorSE(:,n) = sqrt( (pctCor(:,n).*(1-pctCor(:,n))) ./ N(:,n) );
    
    figure(n); set(gcf,'Color',[1 1 1],'Position',[400+10*n 700-10*n 900 300],'PaperPositionMode','auto'); clf;
    subplot(1,3,1);
    errorbar(uCoh,RTmean(:,n),RTse(:,n),'bo-')
        ylabel('Mean RT (ms)')
        xlabel('Motion Coherence (%)')
        title ('RT')

    subplot(1,3,2);
    errorbar(uCoh,pctCor(:,n),pctCorSE(:,n),'ro-')
        ylabel('Percent correct')
        xlabel('Motion Coherence (%)')
        title ('accuracy')

    if ~all(isnan(confMean(:,n)))
        subplot(1,3,3);
        K = ~isnan(confMean(:,n)); % some cohs were not used in conf experiment and thus gave NaN in above loop; omit those
        errorbar(uCoh(K),confMean(K,n),confSE(K,n),'go-')
            ylabel('Mean confidence report')
            xlabel('Motion Coherence (%)')
            title ('confidence')
    end
        
    pause;
 
end
      
  
  
%% for Age, just add a logical statement to the 'I' index variable:

% preallocate, now with the second dimension being age (two possible values)
pctCor = nan(length(uCoh),2);
RTmean = nan(length(uCoh),2);
RTse = nan(length(uCoh),2);
confMean = nan(length(uCoh),2);
confSE = nan(length(uCoh),2);
N = nan(length(uCoh),2);

for a = 1:2 % age index, 1 = >50, 2 = <50
    
    for c = 1:length(uCoh)
        I = goodTrials & data.Coherence==uCoh(c) & data.AGELT50==a;
        N(c,a) = sum(I);
        
        pctCor(c,a) = sum(data.Correct(I)==2) / N(c,a); % 1=err, 2=corr

        RTmean(c,a) = mean(data.ReactionTime(I));
        RTse(c,a) = std(data.ReactionTime(I))/sqrt(N(c,a));
                                                          
        J = I & strcmp(data.Test,'test2');
        confMean(c,a) = nanmean(data.Confidence(J));
        confSE(c,a) = nanstd(data.Confidence(J))/sqrt(sum(J));
    end
    
    pctCorSE(:,a) = sqrt( (pctCor(:,a).*(1-pctCor(:,a))) ./ N(:,a) );
    
end

figure(1000); set(gcf,'Color',[1 1 1],'Position',[500 500 900 300],'PaperPositionMode','auto'); clf;
subplot(1,3,1);
errorbar(uCoh,RTmean(:,1),RTse(:,1),'bo-','MarkerFaceColor','b'); hold on;
errorbar(uCoh,RTmean(:,2),RTse(:,2),'bs--');
    ylabel('Mean RT (ms)')
    xlabel('Motion Coherence (%)')
    title ('RT')
legend('age>50','age<50');
    
subplot(1,3,2);
errorbar(uCoh,pctCor(:,1),pctCorSE(:,1),'ro-','MarkerFaceColor','r'); hold on;
errorbar(uCoh,pctCor(:,2),pctCorSE(:,2),'rs--')
    ylabel('Percent correct')
    xlabel('Motion Coherence (%)')
    title ('accuracy')

subplot(1,3,3);
K = ~isnan(confMean(:,1)); % some cohs were not used in conf experiment and thus gave NaN in above loop; omit those
errorbar(uCoh(K),confMean(K,1),confSE(K,1),'go-','MarkerFaceColor','g'); hold on;
K = ~isnan(confMean(:,2)); % some cohs were not used in conf experiment and thus gave NaN in above loop; omit those
errorbar(uCoh(K),confMean(K,2),confSE(K,2),'gs--');
    ylabel('Mean confidence report')
    xlabel('Motion Coherence (%)')
    title ('confidence')
    
    

%% Repeat for disorder group

% preallocate, now with the second dimension being age (two possible values)
pctCor = nan(length(uCoh),2);
RTmean = nan(length(uCoh),2);
RTse = nan(length(uCoh),2);
confMean = nan(length(uCoh),2);
confSE = nan(length(uCoh),2);
N = nan(length(uCoh),2);

for a = 1:2 % age index, 1 = >50, 2 = <50
    
    for c = 1:length(uCoh)
        I = goodTrials & data.Coherence==uCoh(c) & data.Groupx2==a;
        N(c,a) = sum(I);
        
        pctCor(c,a) = sum(data.Correct(I)==2) / N(c,a); % 1=err, 2=corr

        RTmean(c,a) = mean(data.ReactionTime(I));
        RTse(c,a) = std(data.ReactionTime(I))/sqrt(N(c,a));
                                                          
        J = I & strcmp(data.Test,'test2');
        confMean(c,a) = nanmean(data.Confidence(J));
        confSE(c,a) = nanstd(data.Confidence(J))/sqrt(sum(J));
    end
    
    pctCorSE(:,a) = sqrt( (pctCor(:,a).*(1-pctCor(:,a))) ./ N(:,a) );
    
end

figure(1000); set(gcf,'Color',[1 1 1],'Position',[500 500 900 300],'PaperPositionMode','auto'); clf;
subplot(1,3,1);
errorbar(uCoh,RTmean(:,1),RTse(:,1),'bo-','MarkerFaceColor','b'); hold on;
errorbar(uCoh,RTmean(:,2),RTse(:,2),'bs--');
    ylabel('Mean RT (ms)')
    xlabel('Motion Coherence (%)')
    title ('RT')
legend('control','ocd');
    
subplot(1,3,2);
errorbar(uCoh,pctCor(:,1),pctCorSE(:,1),'ro-','MarkerFaceColor','r'); hold on;
errorbar(uCoh,pctCor(:,2),pctCorSE(:,2),'rs--')
    ylabel('Percent correct')
    xlabel('Motion Coherence (%)')
    title ('accuracy')

subplot(1,3,3);
K = ~isnan(confMean(:,1)); % some cohs were not used in conf experiment and thus gave NaN in above loop; omit those
errorbar(uCoh(K),confMean(K,1),confSE(K,1),'go-','MarkerFaceColor','g'); hold on;
K = ~isnan(confMean(:,2)); % some cohs were not used in conf experiment and thus gave NaN in above loop; omit those
errorbar(uCoh(K),confMean(K,2),confSE(K,2),'gs--');
    ylabel('Mean confidence report')
    xlabel('Motion Coherence (%)')
    title ('confidence')
    
    
    