function dots3DMP_ConfDelta(data,gfit,cohs,deltas,hdgs,conftask)

% 3 analyses

% 1. compare shifts of choice mu vs conf mu in non-conflict and conflict conditions 
% 2. piecewise comparison of PRight for fixed absolute heading, different conflicts as function of coh
% 3. average confidence in conflict vs no conflict (low headings only)

%% compare shifts
n = 1;
for c=1:length(cohs)
    for d = 1:length(deltas)
        choiceBias(n) = gfit.choice.mu(3,c,d);
        choiceBiasSE(n) = gfit.choice.muSE(3,c,d);
        
        confBias(n) = gfit.conf.mu(3,c,d);
        confBiasSE(n) = gfit.conf.muSE(3,c,d);
        
        RTBias(n) = gfit.RT.mu(3,c,d);
        RTBiasSE(n) = gfit.RT.muSE(3,c,d);
        n = n+1;
    end
end
figure(20); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
[hsym,hxe,hye] = errorbar2(choiceBias, confBias, choiceBiasSE, confBiasSE, 'o', 2);
set(hsym,'MarkerSize',10,'MarkerFaceColor','w','Color','k');
hold on; plot([-2 2],[-2 2],'k--','LineWidth',2); axis square;
xlim([-2 2]); ylim([-2 2]);
set(gca,'Xtick',-2:1:2,'Ytick',-2:1:2);
xlabel('Choice shift (deg)'); ylabel('Confidence shift (deg)');
changeAxesFontSize(gca,20,20); set(gca,'box','off')

[B,Bint,R,rint] = regress(confBias',[ones(size(choiceBias)); choiceBias]');

%% 1. relationship bw conf and weights

hdgInd = 3;

% for each coh, bin trials by conf (try median split first), then calc
% weights separately for low and high conf trials
% this logic doesn't work, because confidence in the combined condition per
% se would not map onto weights even under the hypothesis
% BUT we can do it piecewise:

% high coh, left side of the graph, delta +3,
I = data.modality==3 & ismember(data.coherence,cohs(2)) & data.heading<=hdgs(hdgInd) & data.delta==deltas(end);
% when conf is lower, Pright should be higher (anticorrelated)

if conftask==1
    confLow = data.conf(I) < median(data.conf(I));
    confHigh = data.conf(I) >= median(data.conf(I));
else
    confLow = data.PDW(I) == 0;
    confHigh = data.PDW(I) == 1;
end

PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

% clr = {[255 180 0]./255, [175 238 238]./255, [255 140 255]./255};

try clr = cbrewer('qual','Paired',4);
catch
    clr = [0.6510    0.8078    0.8902
            0.1216    0.4706    0.7059
            0.6980    0.8745    0.5412
            0.2000    0.6275    0.1725];
end

figure(21); set(gcf,'Color',[1 1 1],'Position',[50 20 250 200],'PaperPositionMode','auto'); clf;
h = plot([1 2], [PrightLow PrightHigh], 'ok-', 'MarkerSize',10, 'MarkerFaceColor','k');
set(h,'Color', clr(1,:), 'MarkerFaceColor',clr(1,:));
L{1} = sprintf('d=%d, c=%g, h=%g',deltas(end),cohs(2),hdgs(hdgInd));

xlim([0.5 2.5]); ylim([0 1]);
set(gca,'Xtick',[1 2],'XTickLabel',{'Low','High'},'Ytick',0:0.2:1);
ylabel('P(Right)');
changeAxesFontSize(gca,16,16); set(gca,'box','off');
hold on;
% if exportfigs; export_fig('medianSplit1','-eps'); end 

%% low coh, left side of the graph, delta -3, 
I = data.modality==3 & data.coherence==cohs(1) & data.heading<=hdgs(hdgInd) & data.delta==deltas(1);
% when conf lower, Pright should be higher (anticorrelated)
if conftask==1
    confLow = data.conf(I) < median(data.conf(I));
    confHigh = data.conf(I) >= median(data.conf(I));
else
    confLow = data.PDW(I) == 0;
    confHigh = data.PDW(I) == 1;
end

PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'o-', 'MarkerSize', 10);
set(h,'Color', clr(2,:), 'MarkerFaceColor', clr(2,:));
% if exportfigs; export_fig('medianSplit2','-eps'); end
L{2} = sprintf('d=%d, c=%g, h=%g',deltas(1),cohs(1),hdgs(hdgInd));


%% 2. high coh, right side of the graph, delta -3, 
I = data.modality==3 & ismember(data.coherence,cohs(2:end)) & data.heading>=hdgs(end-hdgInd+1) & data.delta==deltas(1);
% when conf is lower, Pright should be lower (correlated)
if conftask==1
    confLow = data.conf(I) < median(data.conf(I));
    confHigh = data.conf(I) >= median(data.conf(I));
else
    confLow = data.PDW(I) == 0;
    confHigh = data.PDW(I) == 1;
end

PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'ob-', 'MarkerSize',10, 'MarkerFaceColor','b');
set(h,'Color', clr(3,:), 'MarkerFaceColor',clr(3,:));
% if exportfigs; export_fig('medianSplit3','-eps'); end
L{3} = sprintf('d=%d, c=%g, h=%g',deltas(1),cohs(2),hdgs(end-hdgInd+1));


%% low coh, right side of the graph, delta +3, 
I = data.modality==3 & data.coherence==cohs(1) & data.heading>=hdgs(end-hdgInd+1) & data.delta==deltas(end);
% when conf lower, Pright should be lower (correlated)
if conftask==1
    confLow = data.conf(I) < median(data.conf(I));
    confHigh = data.conf(I) >= median(data.conf(I));
else
    confLow = data.PDW(I) == 0;
    confHigh = data.PDW(I) == 1;
end

PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'ob-', 'MarkerSize',10, 'MarkerFaceColor','b');
set(h,'Color', clr(4,:), 'MarkerFaceColor',clr(4,:));
% if exportfigs; export_fig('medianSplit4','-eps'); end
L{4} = sprintf('d=%d, c=%g, h=%g',deltas(end),cohs(1),hdgs(end-hdgInd+1));

legend(L);



%% 3. does conflict affect confidence?

maxhdg = 1.5;
for c = 1:length(cohs)
    I = data.modality==3 & data.coherence==cohs(c) & data.delta==0;
    I = I & abs(data.heading)<=maxhdg;
    
    if conftask==1
        confDzero(c) = mean(data.conf(I));
        confDzeroSE(c) = std(data.conf(I))/sqrt(sum(I));
    else
        confDzero(c) = mean(data.PDW(I));
        confDzeroSE(c) = sqrt( (confDzero(c).*(1-confDzero(c))) ./ sum(I));
    end
    
    rtDzero(c) = mean(data.RT(I));
    rtDzeroSE(c) = std(data.RT(I))/sqrt(sum(I));

    J = data.modality==3 & data.coherence==cohs(c) & data.delta~=0;
    J = J & abs(data.heading)<=maxhdg;
    
    if conftask==1
        confDnonzero(c) = mean(data.conf(J));
        confDnonzeroSE(c) = std(data.conf(J))/sqrt(sum(J));
    else
        confDnonzero(c) = mean(data.PDW(J));
        confDnonzeroSE(c) = sqrt( (confDzero(c).*(1-confDzero(c))) ./ sum(J));
    end
    
    rtDnonzero(c) = mean(data.RT(J));
    rtDnonzeroSE(c) = std(data.RT(J))/sqrt(sum(J));

%     [~,pvalConf(c)] = ttest2(data.conf(I),data.conf(J));
%     [~,pvalRT(c)] = ttest2(data.RT(I),data.RT(J));
end

barx = 1:length(cohs);
bary = [confDzero ; confDnonzero]';
errlow = [confDzeroSE ; confDnonzeroSE]';
errhigh = errlow;

figure(809); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
bar(barx,bary); hold on;
er = errorbar([barx-0.14;barx+0.14]',bary,errlow,errhigh,'k.');    
ylim([0.5 1]);
set(gca,'XTickLabel',[10 50 90],'Ytick',0:0.2:1);
xlabel('Visual coherence (%)'); 
if conftask==1
    ylabel('Mean SEP');
else
    ylabel('Mean P(High Bet)')
end
legend(sprintf('\{Delta}=0', '\{Delta}~=0','interpreter','latex')); legend('boxoff');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
% if exportfigs; export_fig('confDeltaZeroNonzero','-eps'); end

%{
bary = [rtDzero ; rtDnonzero]';
errlow = [rtDzeroSE ; rtDnonzeroSE]';
errhigh = errlow;

figure(810); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
bar(barx,bary); hold on;
er = errorbar([barx-0.14;barx+0.14]',bary,errlow,errhigh,'k.');    
set(gca,'XTickLabel',[10 50 90]);
xlabel('Visual coherence (%)'); ylabel('Mean RT (s)');
% legend('no conflict', 'conflict'); legend('boxoff');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
% if exportfigs; export_fig('RTDeltaZeroNonzero','-eps'); end

%}