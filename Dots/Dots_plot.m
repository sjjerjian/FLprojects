%% all trials

figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 500 450 800],'PaperPositionMode','auto'); clf;

if length(cohs)==11
    ticks = cohs([1 2 4 6 8 10 11]); % which cohs to label on the x axis
else
    ticks = cohs;
end

subplot(3,1,1);
plot(xVals,yVals1,'k-'); hold on;
errorbar(cohs, pRight, pRightSE, 'ko');
set(gca,'xtick',ticks,'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion rightward choices');
changeAxesFontSize(gca,11,11);

subplot(3,1,2);
errorbar(cohs, RTmean, RTse, 'bo-');
set(gca,'xtick',ticks,'tickdir','out');
xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
changeAxesFontSize(gca,11,11);

subplot(3,1,3);
errorbar(cohs, pHigh, pHighSE, 'ro-');
set(gca,'xtick',ticks,'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion high bet');
changeAxesFontSize(gca,11,11);

%% separate choice+RT by high/low bet, and p(high) by corr/err

figure(102);
set(gcf,'Color',[1 1 1],'Position',[500 500 450 800],'PaperPositionMode','auto'); clf;

subplot(3,1,1);
plot(xVals,yVals2,'k-',xVals,yVals3,'k--'); hold on;
errorbar(cohs, pRightHigh, pRightSEhigh, 'ko', 'MarkerFaceColor', 'k');
errorbar(cohs, pRightLow, pRightSElow, 'ko', 'MarkerFaceColor', 'w');
set(gca,'xtick',ticks,'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion rightward choices');
legend('High bet','Low bet','Location','Northwest'); legend('boxoff')
changeAxesFontSize(gca,11,11);

subplot(3,1,2);
errorbar(cohs, RTmeanHigh, RTseHigh, 'bo-', 'MarkerFaceColor', 'b'); hold on;
errorbar(cohs, RTmeanLow, RTseLow, 'bo--', 'MarkerFaceColor', 'w');
set(gca,'xtick',ticks,'tickdir','out');
xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
legend('High bet','Low bet','Location','Northwest'); legend('boxoff')
changeAxesFontSize(gca,11,11);

subplot(3,1,3);
errorbar(cohs, pHighCorr, pHighSEcorr, 'ro-', 'MarkerFaceColor', 'r'); hold on;
errorbar(cohs, pHighErr, pHighSEerr, 'ro--', 'MarkerFaceColor', 'w');
set(gca,'xtick',ticks,'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion high bet');
legend('Corr','Err','Location','Southwest'); legend('boxoff')
changeAxesFontSize(gca,11,11);
