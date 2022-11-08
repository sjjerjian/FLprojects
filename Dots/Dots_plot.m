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
    confAxisLabel = 'Confidence (mean SEP or rating)';
    legEntry{1} = 'high conf'; legEntry{2} = 'low conf';
end
if conftask==2
    conf = parsedData.pHigh;
    confSE = parsedData.pHighSE;
    confCorr = parsedData.pHighCorr;
    confSEcorr = parsedData.pHighSEcorr;
    confErr = parsedData.pHighErr;
    confSEerr = parsedData.pHighSEerr;
    confAxisLabel = 'Proportion high bet';
    legEntry{1} = 'high bet'; legEntry{2} = 'low bet';
%    confAxisLabel = 'proportion high conf rating'; % temp, for doubtConf
%    legEntry{1} = 'high conf'; legEntry{2} = 'low conf'; % temp, for doubtConf
end
    
figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;

% which cohs to label on the x axis:
if length(cohs)==12
    ticks = cohs([1 2 4 9 11 12]); 
elseif length(cohs)==11
    ticks = cohs([1 2 4 6 8 10 11]); 
elseif length(cohs)==10 % doubtconf
    ticks = [cohs(1:5)' 0 cohs(6:10)'];
%     xlabels = {'-0.7','-0.45','-0.25','','','0','','','0.25','0.45','0.7'};
    xlabels = {'-0.7','-0.45','','-0.15','','0','','0.15','','0.45','0.7'};
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

xl = [-1.05*max(cohs) 1.05*max(cohs)];
ylConf = [0.4 1]; % hanzo
ylRT = [0.3 0.7]; % hanzo
% ylConf = [0.2 1]; % genji
% ylRT = [0.3 0.9]; % genji

subplot(nplots,1,1);
if wFit==0; plot(parsedData.xVals,parsedData.yVals1,'k-'); hold on; end
errorbar(cohs, parsedData.pRight, parsedData.pRightSE, 'ko'); hold on;
try
    set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
catch
    set(gca,'xtick',ticks,'tickdir','out');
end
ylim([0 1]); xlim(xl);
xlabel('Motion strength (%coh)');
ylabel('Proportion rightward choices');
changeAxesFontSize(gca,14,14);

if RTtask
    subplot(nplots,1,2);
    errorbar(cohs, parsedData.RTmean, parsedData.RTse, ['bo' line1]);
    try
        set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
    catch
        set(gca,'xtick',ticks,'tickdir','out');
    end
    xlim(xl);
    xlabel('Motion strength (%coh)'); ylabel('Reaction time (s)');
    changeAxesFontSize(gca,14,14);
    if conftask
        subplot(nplots,1,3);
        errorbar(cohs, conf, confSE, ['ro' line1]);
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        changeAxesFontSize(gca,14,14);
    end    
else 
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, conf, confSE, ['ro' line1]);
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        changeAxesFontSize(gca,14,14);
    end
end


%% separate RT by corr/err [maybe pick this up later]

% % if RTtask
% %    
% %     figure(102);
% %     set(gcf,'Color',[1 1 1],'Position',[600 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;
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
try
    set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
catch
    set(gca,'xtick',ticks,'tickdir','out');
end
ylim([0 1]); xlim(xl);
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');
legend(legEntry{1},legEntry{2},'Location','Northwest'); legend('boxoff')
changeAxesFontSize(gca,14,14);

if RTtask
    subplot(nplots,1,2);
    errorbar(cohs, parsedData.RTmeanHigh, parsedData.RTseHigh, ['bo' line1], 'MarkerFaceColor', 'b'); hold on;
    errorbar(cohs, parsedData.RTmeanLow, parsedData.RTseLow, ['bo' line2], 'MarkerFaceColor', 'w');
    try
        set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
    catch
        set(gca,'xtick',ticks,'tickdir','out');
    end
    xlim(xl);
    xlabel('Motion strength (%coh)'); ylabel('Reaction time (s)');
    legend(legEntry{1},legEntry{2},'Location','Northwest'); legend('boxoff')
    changeAxesFontSize(gca,14,14);
    if conftask
        subplot(nplots,1,3);
        errorbar(cohs, confCorr, confSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, confErr, confSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        legend('Corr','Err','Location','Southwest'); legend('boxoff')
        changeAxesFontSize(gca,14,14);
    end
else
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, confCorr, confSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, confErr, confSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        try
            set(gca,'xtick',ticks,'tickdir','out','xticklabel',xlabels);
        catch
            set(gca,'xtick',ticks,'tickdir','out');
        end
        ylim([0 1]); xlim(xl);
        xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
        legend('Corr','Err','Location','North'); legend('boxoff')
        changeAxesFontSize(gca,14,14);
    end
end

end


%% nicer figs + indiv panels for talk

if forTalk

% choice
figure(105); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
if wFit==0; plot(parsedData.xVals,parsedData.yVals1,'k-','LineWidth',2); hold on; end
errorbar(cohs, parsedData.pRight, parsedData.pRightSE, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0 1]); xlim(xl);
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
export_fig('dots_choice','-eps');
    % choice split
figure(106); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
if wFit==0
    plot(parsedData.xVals,parsedData.yVals2,'k-','LineWidth',2); hold on;
    plot(parsedData.xVals,parsedData.yVals3,'k--','LineWidth',2);
end
% errorbar(cohs, parsedData.pRightHigh, parsedData.pRightSEhigh, ['ko' line1], 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'LineWidth', 2);
errorbar(cohs, parsedData.pRightHigh, parsedData.pRightSEhigh, 'ko', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
errorbar(cohs, parsedData.pRightLow, parsedData.pRightSElow, 'ko', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
h = legend(legEntry{1},legEntry{2},'Location','northwest');
set(h,'FontSize',16); legend('boxoff');
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0 1]); xlim(xl);
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
if wFit==0; export_fig('dots_choice_split','-eps'); end

if RTtask    
    % RT
    figure(107); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, parsedData.RTmean, parsedData.RTse, ['o' line1], 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    ylim(ylRT); xlim(xl);
    set(gca,'xtick',-0.5:0.25:0.5,'ytick',ylRT(1):0.1:ylRT(2),'tickdir','out','box','off');
    xlabel('Motion strength (%coh)'); ylabel('Reaction time (s)');
    changeAxesFontSize(gca,20,20);
    if wFit==0; export_fig('dots_RT','-eps'); end 
        % RT split
    figure(108); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, parsedData.RTmeanHigh, parsedData.RTseHigh, ['bo' line1], 'MarkerFaceColor', 'b', 'MarkerSize', 10, 'LineWidth', 2); hold on;
    errorbar(cohs, parsedData.RTmeanLow, parsedData.RTseLow, ['bo' line2], 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    set(gca,'xtick',-0.5:0.25:0.5,'ytick',ylRT(1):0.1:ylRT(2),'tickdir','out','box','off');
    xlim(xl); ylim(ylRT);
    xlabel('Motion strength (%coh)'); ylabel('Reaction time (s)');
    changeAxesFontSize(gca,20,20);
%     legend(legEntry{1},legEntry{2},'Location','Northeast'); legend('boxoff')
    if wFit==0; export_fig('dots_RT_split','-eps'); end
end

if conftask
    % PDW
    figure(109); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, conf, confSE, ['o' line1], 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
    ylim(ylConf); xlim(xl);
    xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
    changeAxesFontSize(gca,20,20);
    if wFit==0; export_fig('dots_conf','-eps'); end
        % PDW split
    figure(110); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
    errorbar(cohs, confCorr, confSEcorr, ['o' line1], 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'LineWidth', 2); hold on;
    errorbar(cohs, confErr, confSEerr, ['o' line2], 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
    set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
    ylim(ylConf); xlim(xl);
    xlabel('Motion strength (%coh)'); ylabel(confAxisLabel);
    h = legend('Correct','Error','Location','southeast');
    set(h,'FontSize',16); legend('boxoff');
    changeAxesFontSize(gca,20,20);
    if wFit==0; export_fig('dots_conf_split','-eps'); end
end


end

