

% %% temp: mapping figs for talk
% 
% % open('/Users/chris/Downloads/Chris_LIP_Mapping_1.fig')
% for n=1:3
%     subplot(3,1,n);
%     xlim([10 40]);
%     set(gca,'xtick',10:5:40,'xticklabel',{'-15','-10','-5','0','5','10','15'});
%     set(gca,'ytick',0:5:30,'yticklabel',{'-15','-10','-5','0','5','10','15'});
%     changeAxesFontSize(gca,14,14);
%     title([]);
% end
% 


%% 

clear
close all
load('/Users/chris/Documents/MATLAB/Fall2021_Dataset_LIP.mat');
dataCell = dataLIPCell; clear dataLIPCell;

dbstop if error


% Population vector through time:

% each cell has a center of mass of RF = theta in polar coords
% each cell has a firing rate R at time t
% this corresponds to a vector with amplitude R and angle theta
% add up the vectors of all the cells, and there is a population vector for intended saccade direction (internal state, or decision variable) at every time step
% The prediction is that this should start out close to 0/180, and eventually end at +/- 45 for high/low bet
% 
% so the first thing to do is see if you can do this on individual trials. depends on the sampling of upper and lower field in a given session. Let?s find that out first.


% Step 1: get polar angle of RF center-of-mass for each unit in a given session


% 1 =LeftLow, 2=LeftHigh, 3=RightLow, and 4=RightHigh


% anyway, let's assume we can get that later. Then what.


% Step 2: get FR through time

convKernel = fspecial('average', [1 40]); % N ms wide boxcar (acausal, centered)
minRT = 400; 
maxCoh = 0.128;

% indices for the preexisting PSTH variable (dots onset is at t=300)
tminOn = 200;
tmaxOn = 1000; 
tDotsOn = 300;
tAxisOn = tminOn-tDotsOn:tmaxOn-tDotsOn;

tminOff = 2600;
tmaxOff = 3300; 
tDotsOff = 3200; % dots ofset is at t-minus 300
tAxisOff = tminOff-tDotsOff:tmaxOff-tDotsOff;

contra = 1; % 1=left, 2=right
baseline = 1:tDotsOn+50; % for baseline subtraction

%% 
clear psthOn psthOff psthNormOn psthNormOff RFind RFtheta sessID
C=1;

for n = 1:length(dataCell)

    nU = length(dataCell{n}.Exp.spikeTimes);
    mStart = dataCell{n}.Exp.openEvents.motionStart';
    nTrials = length(mStart);
    mEnd = dataCell{n}.Exp.openEvents.motionEnd';
    RT = round((mEnd-mStart)*1000);
    validTr = RT<=1400; % exclude very long RTs, to rule out 'lapses'

    J = validTr & abs(dataCell{n}.Exp.openEvents.coherence')<=maxCoh & RT >= minRT;
    for c = 1:nU
        
        %bline = mean(nanmean(dataCell{n}.Exp.MotionOn_PSTH{c}(:,baseline))*1e3);
        % or
        bline = 0;
        
        I = J & dataCell{n}.Exp.openEvents.choice'==contra & dataCell{n}.Exp.openEvents.pdw'==0; % contra-low       
        psthOn{1}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOn_PSTH{c}(I,tminOn:tmaxOn))*1e3-bline, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice'==contra & dataCell{n}.Exp.openEvents.pdw'==1; % contra-high
        psthOn{2}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOn_PSTH{c}(I,tminOn:tmaxOn))*1e3-bline, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice'~=contra & dataCell{n}.Exp.openEvents.pdw'==0; % ipsi-low
        psthOn{3}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOn_PSTH{c}(I,tminOn:tmaxOn))*1e3-bline, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice'~=contra & dataCell{n}.Exp.openEvents.pdw'==1; % ipsi-high
        psthOn{4}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOn_PSTH{c}(I,tminOn:tmaxOn))*1e3-bline, convKernel);

        rmax = max([psthOn{1}(C,:) psthOn{2}(C,:) psthOn{3}(C,:) psthOn{4}(C,:)]);
        psthNormOn{1}(C,:) = psthOn{1}(C,:)/rmax;
        psthNormOn{2}(C,:) = psthOn{2}(C,:)/rmax;
        psthNormOn{3}(C,:) = psthOn{3}(C,:)/rmax;
        psthNormOn{4}(C,:) = psthOn{4}(C,:)/rmax;
        
        % repeat for aligned-RT
        I = J & dataCell{n}.Exp.openEvents.choice'==contra & dataCell{n}.Exp.openEvents.pdw'==0; % contra-low      
        psthOff{1}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOff_PSTH{c}(I,tminOff:tmaxOff))*1e3-bline, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice'==contra & dataCell{n}.Exp.openEvents.pdw'==1; % contra-high
        psthOff{2}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOff_PSTH{c}(I,tminOff:tmaxOff))*1e3-bline, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice'~=contra & dataCell{n}.Exp.openEvents.pdw'==0; % ipsi-low
        psthOff{3}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOff_PSTH{c}(I,tminOff:tmaxOff))*1e3-bline, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice'~=contra & dataCell{n}.Exp.openEvents.pdw'==1; % ipsi-high
        psthOff{4}(C,:) = smoothRaster(nanmean(dataCell{n}.Exp.MotionOff_PSTH{c}(I,tminOff:tmaxOff))*1e3-bline, convKernel);
        rmax = max([psthOff{1}(C,:) psthOff{2}(C,:) psthOff{3}(C,:) psthOff{4}(C,:)]);
        psthNormOff{1}(C,:) = psthOff{1}(C,:)/rmax;
        psthNormOff{2}(C,:) = psthOff{2}(C,:)/rmax;
        psthNormOff{3}(C,:) = psthOff{3}(C,:)/rmax;
        psthNormOff{4}(C,:) = psthOff{4}(C,:)/rmax;

        % until we have actual centers of mass, just assign a default theta based
        % on preferred target
        thetas = [130 ; 230];
        RFind(C,1) = dataCell{n}.Mapping.unitType(c);
        RFtheta(C,1) = thetas(RFind(C))*pi/180;
        sessID(C,1) = n;
        C = C+1;
    end
    
end



%% each group of cells in a plot by itself

% high cells, choice=contra, split by wager
figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
plot(tAxisOn,nanmean(psthNormOn{2}(RFind==2,:)),'b-','LineWidth',3); hold on;
plot(tAxisOn,nanmean(psthNormOn{1}(RFind==2,:)),'b--','LineWidth',3);
set(gca, 'XLim', [tAxisOn(1) tAxisOn(end)],'ylim',[0.2 0.7],'ytick',0.2:0.1:0.7,'TickDir','out'); box off;
xlabel('Time from dots on (ms)');
ylabel('Normalized firing rate');
legend('High cells | high bet','High cells | low bet','Location','Northwest');
legend('boxoff');
changeAxesFontSize(gca, 22, 22);
export_fig('HIGH_chContra_wagerSplit', '-eps');

% low cells, choice=contra, split by wager
figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
plot(tAxisOn,nanmean(psthNormOn{2}(RFind==1,:)),'r-','LineWidth',3); hold on;
plot(tAxisOn,nanmean(psthNormOn{1}(RFind==1,:)),'r--','LineWidth',3);
set(gca, 'XLim', [tAxisOn(1) tAxisOn(end)],'ylim',[0.2 0.7],'ytick',0.2:0.1:0.7,'TickDir','out'); box off;
xlabel('Time from dots on (ms)');
ylabel('Normalized firing rate');
legend('Low cells | high bet','Low cells | low bet','Location','Northwest');
legend('boxoff');
changeAxesFontSize(gca, 22, 22);
export_fig('LOW_chContra_wagerSplit', '-eps');


% % 
% % % TEMP
% % 
% % % all cells, choice=contra, split by wager
% % figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
% % plot(tAxisOn,nanmean(psthNormOn{2}(:,:)),'b-','LineWidth',2); hold on;
% % plot(tAxisOn,nanmean(psthNormOn{1}(:,:)),'b--','LineWidth',2);
% % set(gca, 'XLim', [tAxisOn(1) tAxisOn(end)],'ylim',[0.2 0.7],'ytick',0.2:0.1:0.7,'TickDir','out'); box off;
% % xlabel('Time from dots on (ms)');
% % ylabel('Normalized firing rate');
% % legend('All cells | high bet','All cells | low bet','Location','Northwest');
% % legend('boxoff');
% % changeAxesFontSize(gca, 22, 22);
% % export_fig('ALL_chContra_wagerSplit', '-eps');
% % 



% high cells, bet=high, split by choice
figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
plot(tAxisOn,nanmean(psthNormOn{2}(RFind==2,:)),'c-','LineWidth',3); hold on;
plot(tAxisOn,nanmean(psthNormOn{4}(RFind==2,:)),'c--','LineWidth',3);
set(gca, 'XLim', [tAxisOn(1) tAxisOn(end)],'ylim',[0.2 0.7],'ytick',0.2:0.1:0.7,'TickDir','out'); box off;
xlabel('Time from dots on (ms)');
ylabel('Normalized firing rate');
legend('High cells | Contra choice','High cells | Ipsi choice','Location','Northwest');
legend('boxoff');
changeAxesFontSize(gca, 22, 22);
export_fig('HIGH_wagerHigh_chSplit', '-eps');

% low cells, bet=low, split by choice
figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
plot(tAxisOn,nanmean(psthNormOn{1}(RFind==1,:)),'m-','LineWidth',3); hold on;
plot(tAxisOn,nanmean(psthNormOn{3}(RFind==1,:)),'m--','LineWidth',3);
set(gca, 'XLim', [tAxisOn(1) tAxisOn(end)],'ylim',[0.2 0.7],'ytick',0.2:0.1:0.7,'TickDir','out'); box off;
xlabel('Time from dots on (ms)');
ylabel('Normalized firing rate');
legend('Low cells | Contra choice','Low cells | Ipsi choice','Location','Northwest');
legend('boxoff');
changeAxesFontSize(gca, 22, 22);
export_fig('LOW_wagerHigh_chSplit', '-eps');




%% repeat for aligned-RT
% 
% figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
% plot(tAxisOff,nanmean(psthNormOff{1}(RFind==1,:)),'b-','LineWidth',2); hold on;
% plot(tAxisOff,nanmean(psthNormOff{2}(RFind==1,:)),'b--','LineWidth',2);
% plot(tAxisOff,nanmean(psthNormOff{1}(RFind==2,:)),'r-','LineWidth',2);
% plot(tAxisOff,nanmean(psthNormOff{2}(RFind==2,:)),'r--','LineWidth',2);
% yl = ylim;
% plot([0 0],yl,'k--');
% set(gca, 'XLim', [tAxisOff(1) tAxisOff(end)],'YAxisLocation', 'Right','TickDir','out'); box off;
% xlabel('Time from saccade (ms)');
% ylabel('Normalized firing rate');
% legend('High cells | high bet','High cells | low bet','Low cells | high bet','Low cells | low bet','Location','Northwest');
% legend('boxoff');
% changeAxesFontSize(gca, 22, 22);



% %%
% for p = 1:2
%     for t = 1:Tmax
%         [X,Y] = pol2cart(RFtheta,psth{p}(:,t));
%         [vecsum_theta(p,t), vecsum_R(p,t)] = cart2pol(sum(X),sum(Y));
%     end
%     figure(p);
%     subplot(2,1,1); plot(vecsum_theta(p,:)*180/pi); ylabel('theta');
%     subplot(2,1,2); plot(vecsum_R(p,:)); ylabel('R');
% 
%     figure(p*10);
%     for t = 400:5:Tmax
%         startTime = max([400 t-10]);
%         polarplot(vecsum_theta(p,startTime:t)*180/pi+360,vecsum_R(p,startTime:t),'b.-');
%         rlim([0 180]);
%         title(['t=' num2str(t)]);
%         pause(0.2); drawnow;        
%     end
% end




% figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
% plot(tAxisOn,nanmean(psthNormOn{1}(RFind==1,:)),'b-','LineWidth',2); hold on;
% plot(tAxisOn,nanmean(psthNormOn{2}(RFind==2,:)),'r-','LineWidth',2);
% set(gca, 'XLim', [tAxisOn(1) tAxisOn(end)], 'TickDir','out'); box off;
% xlabel('Time from dots on (ms)');
% ylabel('Normalized firing rate');
% legend('High cells | high bet','Low cells | low bet','Location','Northwest');
% legend('boxoff');
% changeAxesFontSize(gca, 22, 22);
% 
% figure; set(gcf,'Color',[1 1 1],'Position',[20 440 800 600],'PaperPositionMode','auto');
% plot(tAxisOn,nanmean(psthNormOn{2}(RFind==1,:)),'b-','LineWidth',2); hold on;
% plot(tAxisOn,nanmean(psthNormOn{1}(RFind==2,:)),'r-','LineWidth',2);
% set(gca, 'XLim', [tAxisOn(1) tAxisOn(end)], 'TickDir','out'); box off;
% xlabel('Time from dots on (ms)');
% ylabel('Normalized firing rate');
% legend('High cells | low bet','Low cells | high bet','Location','Northwest');
% legend('boxoff');
% changeAxesFontSize(gca, 22, 22);


