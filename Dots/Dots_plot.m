function Dots_plot(parsedData,cohs,conftask,RTtask,wFit,forTalk)


if nargin<5, wFit=0; end
if nargin<6, forTalk=0; end

nplots = 1+RTtask+double(conftask>0);

if conftask==1
    conf = parsedData.confMean;
    confSE = parsedData.confSE;
    confCorr = parsedData.confMeanCorr;
    confSEcorr = parsedData.confSEcorr;
    confErr = parsedData.confMeanErr;
    confSEerr = parsedData.confSEerr;
end
if conftask==2
    conf = parsedData.pHigh;
    confSE = parsedData.pHighSE;
    confCorr = parsedData.pHighCorr;
    confSEcorr = parsedData.pHighSEcorr;
    confErr = parsedData.pHighErr;
    confSEerr = parsedData.pHighSEerr;
end
    
figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;

if length(cohs)==11
    ticks = cohs([1 2 4 6 8 10 11]); % which cohs to label on the x axis
else
    ticks = cohs;
end

if wFit
    line1=[];
    line2=[];
else
    line1='-';
    line2='--';
end

subplot(nplots,1,1);
if wFit==0; plot(parsedData.xVals,parsedData.yVals1,'k-'); hold on; end
errorbar(cohs, parsedData.pRight, parsedData.pRightSE, 'ko'); hold on;
set(gca,'xtick',ticks,'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion rightward choices');
changeAxesFontSize(gca,11,11);

if RTtask
    subplot(nplots,1,2);
    errorbar(cohs, parsedData.RTmean, parsedData.RTse, ['bo' line1]);
    set(gca,'xtick',ticks,'tickdir','out');
    xlim([-0.55 0.55]);
    xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
    changeAxesFontSize(gca,11,11);
    if conftask
        subplot(nplots,1,3);
        errorbar(cohs, conf, confSE, ['ro' line1]);
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        changeAxesFontSize(gca,11,11);
    end    
else 
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, conf, confSE, ['ro' line1]);
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        changeAxesFontSize(gca,11,11);
    end
end


%% separate RT by corr/err  [maybe pick this up later]

% % if RTtask
% %    
% %     figure(102);
% %     set(gcf,'Color',[1 1 1],'Position',[600 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;
% % 
% % 
% % end


%% separate choice+RT by high/low bet, and p(high) (or conf rating) by corr/err

if conftask

figure(102);
set(gcf,'Color',[1 1 1],'Position',[600 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;

subplot(nplots,1,1);
if wFit==0; plot(parsedData.xVals,parsedData.yVals2,'k-',parsedData.xVals,parsedData.yVals3,'k--'); hold on; end
errorbar(cohs, parsedData.pRightHigh, parsedData.pRightSEhigh, 'ko', 'MarkerFaceColor', 'k'); hold on;
errorbar(cohs, parsedData.pRightLow, parsedData.pRightSElow, 'ko', 'MarkerFaceColor', 'w');
set(gca,'xtick',ticks,'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion rightward choices');
legend('High bet','Low bet','Location','Northwest'); legend('boxoff')
changeAxesFontSize(gca,11,11);

if RTtask
    subplot(nplots,1,2);
    errorbar(cohs, parsedData.RTmeanHigh, parsedData.RTseHigh, ['bo' line1], 'MarkerFaceColor', 'b'); hold on;
    errorbar(cohs, parsedData.RTmeanLow, parsedData.RTseLow, ['bo' line2], 'MarkerFaceColor', 'w');
    set(gca,'xtick',ticks,'tickdir','out');
    xlim([-0.55 0.55]);
    xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
    legend('High bet','Low bet','Location','Northwest'); legend('boxoff')
    changeAxesFontSize(gca,11,11);
    if conftask
        subplot(nplots,1,3);
        errorbar(cohs, confCorr, confSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, confErr, confSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        legend('Corr','Err','Location','Southwest'); legend('boxoff')
        changeAxesFontSize(gca,11,11);
    end
else
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, confCorr, confSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, confErr, confSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        legend('Corr','Err','Location','Southwest'); legend('boxoff')
        changeAxesFontSize(gca,11,11);
    end
end

end


%% nicer figs + indiv panels for talk

if forTalk

% choice
figure(105); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
plot(parsedData.xVals,parsedData.yVals1,'k-','LineWidth',2); hold on;
errorbar(cohs, parsedData.pRight, parsedData.pRightSE, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
export_fig('dots_pmf','-eps');
    % choice split
figure(106); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
plot(parsedData.xVals,parsedData.yVals2,'r-','LineWidth',2); hold on;
plot(parsedData.xVals,parsedData.yVals3,'b-','LineWidth',2);
errorbar(cohs, parsedData.pRightHigh, parsedData.pRightSEhigh, 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'LineWidth', 2);
errorbar(cohs, parsedData.pRightLow, parsedData.pRightSElow, 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 10, 'LineWidth', 2);
h = legend('High bet','Low bet','Location','northwest');
set(h,'FontSize',16); legend('boxoff');
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
export_fig('dots_pmf_split','-eps');

if RTtask    
    % RT
    figure(107); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, parsedData.RTmean, parsedData.RTse, 'o-', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    ylim([0.3 0.7]); xlim([-0.55 0.55]);
    set(gca,'xtick',-0.5:0.25:0.5,'ytick',0.3:0.1:0.7,'tickdir','out','box','off');
    xlabel('Motion strength (coh)'); ylabel('Reaction time (s)');
    changeAxesFontSize(gca,20,20);
    export_fig('dots_RT','-eps');
        % RT split
    figure(108); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, parsedData.RTmeanHigh, parsedData.RTseHigh, 'ro-', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'LineWidth', 2); hold on;
    errorbar(cohs, parsedData.RTmeanLow, parsedData.RTseLow, 'bo-', 'MarkerFaceColor', 'b', 'MarkerSize', 10, 'LineWidth', 2);
    set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
    xlim([-0.55 0.55]);
    xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
    changeAxesFontSize(gca,20,20);
    legend('High bet','Low bet','Location','Northeast'); legend('boxoff')
    export_fig('dots_RT_split','-eps');
end

if confTask
    % PDW
    figure(109); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, conf, confSE, 'o-', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
    ylim([0.5 1]); xlim([-0.55 0.55]);
    xlabel('Motion strength (coh)'); ylabel('Proportion high bet');
    changeAxesFontSize(gca,20,20);
    export_fig('dots_conf','-eps');
        % PDW split
    figure(110); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, confCorr, confSEcorr, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2); hold on;
% %     errorbar(cohs(2:end-1), confErr(2:end-1), confSEerr(2:end-1), 'o--', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    errorbar(cohs, confErr, confSEerr, 'o--', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
    ylim([0.4 1]); xlim([-0.55 0.55]);
    xlabel('Motion strength (coh)'); ylabel('Proportion high bet');
    h = legend('Correct','Error','Location','southeast');
    set(h,'FontSize',16); legend('boxoff');
    set(h,'Position',[0.6595    0.8328    0.2444    0.1203]);
    changeAxesFontSize(gca,20,20);
    export_fig('dots_conf_split','-eps');
end


end

