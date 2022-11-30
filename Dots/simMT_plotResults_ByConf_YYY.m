% Plot Choice, RT, PDW distributions in Black and white
% Choice logistics are fit seperately to Lo vs Hi confidence trials
% TO DO: Split RT by Hi/Lo confidence

pRightLo = nan(length(cohs),1);
pRightHi = nan(length(cohs),1);
pHigh = nan(length(cohs),1);

% RTmeanHi = nan(length(cohs),1);
% RTmeanLo = nan(length(cohs),1);
% RTLo_SE = nan(length(cohs),1);
% RTHi_SE = nan(length(cohs),1);

RTmean = nan(length(cohs),1);
RT_SE = nan(length(cohs),1);
nc = nan(length(cohs),1);
n2 = nc;
n3 = nc;

for c = 1:length(cohs)

    I = coh(1:size(Ry,1))==cohs(c);   % I = logical index to trials of coh 'c'

    pRLo_J = I & pdw(1:size(Ry,1))==0;
    n3(c) = sum(pRLo_J); % How many times, Lo on coh 'c' trials
    pRightLo(c,1) = sum(pRLo_J & choice(1:size(Ry,1))==1) /  n3(c);
    
    pRHi_J = I & pdw(1:size(Ry,1))==1;
    n2(c) = sum(pRHi_J); % How many times, Hi on coh 'c' trials
    pRightHi(c,1) = sum(pRHi_J & choice(1:size(Ry,1))==1) /  n2(c);
    
    nc(c) = sum(I);
    pHigh(c,1) = sum(pRHi_J) / nc(c); % Same but for Lo choice trials

    RTmean(c) = mean(RT(I));
    RT_SE(c) = std(RT(I)) / sqrt(sum(I));

%     Lo_J = I & pdw(1:size(Ry,1))==0;
%     RTmeanLo(c) = mean(RT(Lo_J)); % mean RT of 'c' trials with Lo bets
% 
%     Hi_J = I & pdw(1:size(Ry,1))==1;
%     RTmeanHi(c) = mean(RT(Hi_J)); % mean RT of 'c' trials with Hi bets
% 
%     RTLo_SE(c) = std(RT(Lo_J)) / sqrt(sum(Lo_J)); % Std error of RT on 'c' trials
%     RTHi_SE(c) = std(RT(Hi_J)) / sqrt(sum(Hi_J)); % Std error of RT on 'c' trials
end

% Standard Error values, for errorbars
pRLo_SE = sqrt(pRightLo.*(1-pRightLo)./n3); % Std error of RLo choice prob at each coh level
pRHi_SE = sqrt(pRightHi.*(1-pRightHi)./n2); % Std error of RHi choice prob at each coh level
pHighSE = sqrt(pHigh.*(1-pHigh)./nc); % For P(Hi bet)

xVals = min(coh):0.01:max(coh); % Coherence values on xAxis

% Plotting Choice
% fit logistic to HIGH bet only
I = pdw==1;
X = coh(I);
y = choice(I)==1;
[B2, ~, stats2] = glmfit(X, y, 'binomial');
yVals2 = glmval(B2,xVals,'logit');
% fit logistic to LOW bet only
I = pdw==0;
X = coh(I);
y = choice(I)==1;
[B3, ~, stats3] = glmfit(X, y, 'binomial');
yVals3 = glmval(B3,xVals,'logit');
% Plot Hi vs Lo confidence choices
figure(10)
set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
plot(xVals,yVals2,'k-','LineWidth',1.5); hold on;
plot(xVals,yVals3,'k--','LineWidth',1.5);
errorbar(cohs, pRightHi, pRHi_SE, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 9, 'LineWidth', 1);
errorbar(cohs, pRightLo, pRLo_SE, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1);
h = legend('High bet','Low bet','Location','northwest');
set(h,'FontSize',16); legend('boxoff');
set(gca,'xtick',-0.5:0.25:0.5,'tickdir','out','box','off');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20); hold on;

% Plotting RT 
figure(20)
errorbar(cohs, RTmean/1000, RT_SE/1000, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 9, 'LineWidth', 1.5);
% errorbar(cohs, RTmeanHi/1000, RTHi_SE/1000, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 9, 'LineWidth', 1.5);hold on;
% errorbar(cohs, RTmeanLo/1000, RTLo_SE/1000, 'o--', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.5);
% ylim([0.3 0.7]); xlim([-0.55 0.55]);
set(gca,'xtick',-0.5:0.25:0.5,'ytick',0.3:0.1:1.25,'tickdir','out','box','off');
xlabel('Motion strength (coh)'); ylabel('Reaction time (s)');
changeAxesFontSize(gca,20,20);

% Plotting Pdw
figure(30)
errorbar(cohs, pHigh, pHighSE, 'o-', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 9, 'LineWidth', 1.5);
set(gca,'xtick',-0.5:0.25:0.5,'ytick',.4:.1:1, 'tickdir','out','box','off');
ylim([0.5 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (coh)'); ylabel('Proportion high bet');
changeAxesFontSize(gca,20,20);





