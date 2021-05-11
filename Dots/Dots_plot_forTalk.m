%% choice

% Dots_fit_cgauss % fix this later

useGauss = 0;

figure(105); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
plot(xVals,yVals1,'k-','LineWidth',2); hold on;
errorbar(cohs, pRight, pRightSE, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
export_fig('dots_pmf','-eps');

figure(106); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
plot(xVals,yVals2,'k-','LineWidth',2); hold on;
plot(xVals,yVals3,'k--','LineWidth',2);
errorbar(cohs, pRightHigh, pRightSEhigh, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
errorbar(cohs, pRightLow, pRightSElow, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
h = legend('High bet','Low bet','Location','northwest');
set(h,'FontSize',16); legend('boxoff');
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
export_fig('dots_pmf_split','-eps');

%% RT

figure(107); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
if useGauss
    beta = [amplRT muRT sigmaRT baselineRT];
    h = plot(xVals, gauss(beta,xVals),'b-','LineWidth',2); hold on;
    errorbar(cohs, RTmean, RTse, 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
else
    errorbar(cohs, RTmean, RTse, 'o-', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
end
ylim([0.3 0.7]); xlim([-0.55 0.55]);
set(gca,'xtick',-0.5:0.25:0.5,'ytick',0.3:0.1:0.7,'tickdir','out','box','off');
xlabel('Motion strength (coh)'); ylabel('Reaction time (s)');
changeAxesFontSize(gca,20,20);
export_fig('dots_RT','-eps');


% split by corr/err or high/low? or both? meh, neither, for now

% % % figure(108); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
% % % errorbar(cohs, pHighCorr, pHighSEcorr, 'o-', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'LineWidth', 2); hold on;
% % % errorbar(cohs, pHighErr, pHighSEerr, 'o--', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
% % % set(gca,'xtick',-0.5:0.25:0.50,'tickdir','out','box','off');
% % % ylim([0.4 1]); xlim([-0.55 0.55]);
% % % xlabel('Motion strength (% coh)'); ylabel('Proportion high bet');
% % % h = legend('Correct','Error','Location','southeast');
% % % set(h,'FontSize',16); legend('boxoff');
% % % set(h,'Position',[0.6595    0.8328    0.2444    0.1203]);
% % % changeAxesFontSize(gca,20,20);
% % % export_fig('dots_RT_split','-eps');
% % % 
% % % subplot(3,1,2);
% % % errorbar(cohs, RTmeanHigh, RTseHigh, 'bo-', 'MarkerFaceColor', 'b'); hold on;
% % % errorbar(cohs, RTmeanLow, RTseLow, 'bo--', 'MarkerFaceColor', 'w');
% % % set(gca,'xtick',ticks,'tickdir','out');
% % % xlim([-0.55 0.55]);
% % % xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
% % % legend('High bet','Low bet','Location','Northwest'); legend('boxoff')
% % % changeAxesFontSize(gca,11,11);



%% PDW

figure(109); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
if useGauss
    beta = [amplConf muConf sigmaConf baselineConf];
    h = plot(xVals, flippedGauss(beta,xVals),'r-','LineWidth',2); hold on;
    errorbar(cohs, pHigh, pHighSE, 'o', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
else
    errorbar(cohs, pHigh, pHighSE, 'o-', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
end
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0.5 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion high bet');
changeAxesFontSize(gca,20,20);
export_fig('dots_conf','-eps');

figure(110); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
errorbar(cohs, pHighCorr, pHighSEcorr, 'o-', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'LineWidth', 2); hold on;
errorbar(cohs, pHighErr, pHighSEerr, 'o--', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0.4 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion high bet');
h = legend('Correct','Error','Location','southeast');
set(h,'FontSize',16); legend('boxoff');
set(h,'Position',[0.6595    0.8328    0.2444    0.1203]);
changeAxesFontSize(gca,20,20);
export_fig('dots_conf_split','-eps');

