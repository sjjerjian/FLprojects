function Dots_plot(parsedData,cohs,conftask,RTtask,wFit)


if nargin<5, wFit=0; end

nplots = 1+RTtask+double(conftask>0);

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
        errorbar(cohs, parsedData.pHigh, parsedData.pHighSE, ['ro' line1]);
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0.5 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        changeAxesFontSize(gca,11,11);
    end    
else 
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, parsedData.pHigh, parsedData.pHighSE, ['ro' line1]);
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0.5 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        changeAxesFontSize(gca,11,11);
    end
end


%% separate choice+RT by high/low bet, and p(high) by corr/err

figure(102);
set(gcf,'Color',[1 1 1],'Position',[500 500 450 200+200*nplots],'PaperPositionMode','auto'); clf;

subplot(nplots,1,1);
if wFit==0 plot(parsedData.xVals,parsedData.yVals2,'k-',parsedData.xVals,parsedData.yVals3,'k--'); hold on; end
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
        errorbar(cohs, parsedData.pHighCorr, parsedData.pHighSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, parsedData.pHighErr, parsedData.pHighSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        legend('Corr','Err','Location','Southwest'); legend('boxoff')
        changeAxesFontSize(gca,11,11);
    end
else
    if conftask
        subplot(nplots,1,2);
        errorbar(cohs, parsedData.pHighCorr, parsedData.pHighSEcorr, ['ro' line1], 'MarkerFaceColor', 'r'); hold on;
        errorbar(cohs, parsedData.pHighErr, parsedData.pHighSEerr, ['ro' line2], 'MarkerFaceColor', 'w');
        set(gca,'xtick',ticks,'tickdir','out');
        ylim([0 1]); xlim([-0.55 0.55]);
        xlabel('motion strength (% coh)'); ylabel('proportion high bet');
        legend('Corr','Err','Location','Southwest'); legend('boxoff')
        changeAxesFontSize(gca,11,11);
    end
end


