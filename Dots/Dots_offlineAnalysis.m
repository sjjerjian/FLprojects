% offline analysis wrapper for PLDAPS data, Dots protocols
% CF spawned it from dots3DMP_offlineAnalysis, begun 10-3-18

clear all; close all

% these will specify which (previously saved) mat file to load
subject = 'hanzo';

% % RT+PDW same session
% dateRange = 20200728:20200730;

% all Hanzo good data
% dateRange = 20190401:20200315;
% OR
% dateRange = 20200101:20200731; % nope nope nope


% Hanzo for R01
% dateRange = 20200728:20200930; % great PDW; misnamed (not really up to 9/30)
% dateRange = 20200728:20200922;
% dateRange = 20200810:20200922;
% dateRange = 20200820:20200922; % really nice RT
dateRange = 20200901:20200922; % good for both, with TrNo<600 or even 800
% dateRange = 20200911:20200922;
% dateRange = 20200917:20200922; 

folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
load([folder file], 'data');

%%

% not doing anyhting with these for now:
try
    data = rmfield(data,{'forgivePDW','trainingSeqPDW','confidenceTime'});
catch
end

% remove invalid trials (fixation breaks, etc)
removethese = isnan(data.choice) | isnan(data.PDW); % may need to add conditions here
% removethese = isnan(data.choice); % not all sessions have PDW! nans will be removed below
% removethese = isnan(data.choice) | ~isnan(data.PDW) % and for RT only, may want to exclude PDW trials
% removethese = isnan(data.choice) | isinf(data.RT) | data.RT==0; % RT, alt (for R01 fig)
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end


% define signed coherence; positive = rightward (0 deg)
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);


% random fixes
data.scoh(abs(data.scoh)<0.001)=0;

% quick look at blocks, for when some need to be excluded
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
data.trialNum = zeros(size(data.filename));
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
    % also number trials within a block, to look at early vs late trials
    data.trialNum(strcmp(blocks(u),data.filename)) = 1:sum(strcmp(blocks(u),data.filename))';
end

% we can be pretty sure blocks with <N trials (say, 10) are to be discarded
removethese = ismember(data.filename,blocks(nTrialsByBlock<10));

% remove early/late trials
removethese = removethese | data.trialNum>600;

fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end



%% parse data

cohs = unique(data.scoh);
dirs = unique(data.direction);

n = nan(length(cohs),1);

n1 = n; % all
n2 = n; % high
n3 = n; % low
nRT1 = n; % all
nRT2 = n; % high
nRT3 = n; % low
nPDW1 = n; % all
nPDW2 = n; % corr
nPDW3 = n; % err

pRight = n1;
pRightHigh = n2;
pRightLow = n3;

RTmean = nRT1;
RTse = nRT1;
RTmeanHigh = nRT2;
RTseHigh = nRT2;
RTmeanLow = nRT3;
RTseLow = nRT3;

pHigh = nPDW1;
pHighCorr = nPDW2;
pHighErr = nPDW3;

% for logistic regression (or whatever)
xVals = linspace(cohs(1),cohs(end),100);

for c = 1:length(cohs)
    % choice
    J = data.scoh==cohs(c);
    n1(c) = sum(J);
    pRight(c) = sum(J & data.choice==2) / n1(c); % 1 is left, 2 is right
    
    JJ = data.scoh==cohs(c) & data.PDW==1;
    n2(c) = sum(JJ);
    pRightHigh(c) = sum(JJ & data.choice==2) / n2(c); 
    
    JJJ = data.scoh==cohs(c) & data.PDW==0;
    n3(c) = sum(JJJ);
    pRightLow(c) = sum(JJJ & data.choice==2) / n3(c); 
    
    % RT
    K = ~isnan(data.RT); % need to start tagging non-RT trials as NaN!
    nRT1(c) = sum(J & K);
    RTmean(c) = mean(data.RT(J & K));
    RTse(c) = std(data.RT(J & K))/sqrt(nRT1(c));
    
    nRT2(c) = sum(JJ & K);
    RTmeanHigh(c) = mean(data.RT(JJ & K));
    RTseHigh(c) = std(data.RT(JJ & K))/sqrt(nRT2(c));

    nRT3(c) = sum(JJJ & K);
    RTmeanLow(c) = mean(data.RT(JJJ & K));
    RTseLow(c) = std(data.RT(JJJ & K))/sqrt(nRT3(c));
    
    % conf
    L = ~isnan(data.PDW);
    nPDW1(c) = sum(J & L);
    pHigh(c) = sum(J & L & data.PDW==1) / nPDW1(c); % 1 is high-bet

    LL = ~isnan(data.PDW) & data.correct==1;
    nPDW2(c) = sum(J & LL);
    pHighCorr(c) = sum(J & LL & data.PDW==1) / nPDW2(c);

    LLL = ~isnan(data.PDW) & data.correct==0;
    nPDW3(c) = sum(J & LLL);
    pHighErr(c) = sum(J & LLL & data.PDW==1) / nPDW3(c);
end

pRightSE = sqrt( (pRight.*(1-pRight)) ./ n1 );
pRightSEhigh = sqrt( (pRightHigh.*(1-pRightHigh)) ./ n2 );
pRightSElow = sqrt( (pRightLow.*(1-pRightLow)) ./ n3 );

pHighSE = sqrt( (pHigh.*(1-pHigh)) ./ nPDW1 );
pHighSEcorr = sqrt( (pHighCorr.*(1-pHighCorr)) ./ nPDW2 );
pHighSEerr = sqrt( (pHighErr.*(1-pHighErr)) ./ nPDW3 );


% fit logistic regression
X = data.scoh;
y = data.choice==2; % 2 is right
[B1, ~, stats1] = glmfit(X, y, 'binomial');

yVals1 = glmval(B1,xVals,'logit');

I = data.PDW==1;
X = data.scoh(I);
y = data.choice(I)==2; % 2 is right
[B2, ~, stats2] = glmfit(X, y, 'binomial');
yVals2 = glmval(B2,xVals,'logit');

I = data.PDW==0;
X = data.scoh(I);
y = data.choice(I)==2; % 2 is right
[B3, ~, stats3] = glmfit(X, y, 'binomial');
yVals3 = glmval(B3,xVals,'logit');

% t test on slopes:
mu = abs(B2(2)-B3(2));
se = sqrt(stats2.se(2)^2 + stats3.se(2)^2);
t = mu/se;
df = sum(~isnan(data.PDW))-length(B2); 
pval_ttest_slopes = 2*(1-tcdf(t,df)) % two-tailed


% t test on pHigh, corr vs err
M = ~isnan(data.PDW) & data.correct==1;
pHighCorr_all = sum(M & data.PDW==1) / sum(M);
MM = ~isnan(data.PDW) & data.correct==0;
pHighErr_all = sum(MM & data.PDW==1) / sum(MM);

pHighSEcorr_all = sqrt( (pHighCorr_all.*(1-pHighCorr_all)) ./ sum(M) );
pHighSEerr_all = sqrt( (pHighErr_all.*(1-pHighErr_all)) ./ sum(MM) );

mu = abs(pHighCorr_all-pHighErr_all);
se = sqrt(pHighSEcorr_all^2 + pHighSEerr_all^2);
t = mu/se; 
df = sum(~isnan(data.PDW))-length(B2); 
pval_ttest_pHighCorrErr = 2*(1-tcdf(t,df)) % two-tailed
% should really use a two-sample Z test for proportions...



%% plot

figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 500 450 800],'PaperPositionMode','auto'); clf;

subplot(3,1,1);
plot(xVals,yVals1,'k-'); hold on;
errorbar(cohs, pRight, pRightSE, 'ko');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion rightward choices');
changeAxesFontSize(gca,14,14);

subplot(3,1,2);
errorbar(cohs, RTmean, RTse, 'bo-');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
changeAxesFontSize(gca,14,14);

subplot(3,1,3);
errorbar(cohs, pHigh, pHighSE, 'ro-');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
ylim([0.5 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion high bet');
changeAxesFontSize(gca,14,14);

% separated by high/low and corr/err
figure(102);
set(gcf,'Color',[1 1 1],'Position',[500 500 450 800],'PaperPositionMode','auto'); clf;

subplot(3,1,1);
plot(xVals,yVals2,'k-',xVals,yVals3,'k--'); hold on;
errorbar(cohs, pRightHigh, pRightSEhigh, 'ko', 'MarkerFaceColor', 'k');
errorbar(cohs, pRightLow, pRightSElow, 'ko', 'MarkerFaceColor', 'w');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion rightward choices');
changeAxesFontSize(gca,14,14);

subplot(3,1,2);
errorbar(cohs, RTmeanHigh, RTseHigh, 'bo-', 'MarkerFaceColor', 'b'); hold on;
errorbar(cohs, RTmeanLow, RTseLow, 'bo--', 'MarkerFaceColor', 'w');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('reaction time (s)');
changeAxesFontSize(gca,14,14);

subplot(3,1,3); % !!WHOAH super x-shape!
errorbar(cohs, pHighCorr, pHighSEcorr, 'ro-', 'MarkerFaceColor', 'r'); hold on;
errorbar(cohs, pHighErr, pHighSEerr, 'ro--', 'MarkerFaceColor', 'w');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion high bet');
changeAxesFontSize(gca,14,14);


%% nicer looking graphs

% Dots_fit_cgauss % fix this later
useGauss = 0;

figure(105); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
plot(xVals*100,yVals1,'k-','LineWidth',2); hold on;
errorbar(cohs*100, pRight, pRightSE, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
set(gca,'xtick',-50:25:50,'tickdir','out','box','off');
ylim([0 1]); xlim([-55 55]);
xlabel('Motion strength (% coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
export_fig('dots_pmf','-eps');

figure(106); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
plot(xVals*100,yVals2,'k-','LineWidth',2); hold on;
plot(xVals*100,yVals3,'k--','LineWidth',2);
errorbar(cohs*100, pRightHigh, pRightSEhigh, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
errorbar(cohs*100, pRightLow, pRightSElow, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
h = legend('High bet','Low bet','Location','northwest');
set(h,'FontSize',14); legend('boxoff');
set(gca,'xtick',-50:25:50,'tickdir','out','box','off');
ylim([0 1]); xlim([-55 55]);
xlabel('Motion strength (% coh)'); ylabel('Proportion rightward choices');
changeAxesFontSize(gca,20,20);
export_fig('dots_pmf_split','-eps');


figure(107); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
if useGauss
    beta = [amplRT muRT sigmaRT baselineRT];
    h = plot(xVals, gauss(beta,xVals),'b-','LineWidth',2); hold on;
    errorbar(cohs, RTmean, RTse, 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
else
    errorbar(cohs, RTmean, RTse, 'o-', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
end
ylim([0.3 0.8]); xlim([-0.55 0.55]);
set(gca,'xtick',-0.50:0.25:0.50,'ytick',0.3:0.1:0.8,'tickdir','out','box','off');
xlabel('Motion strength (% coh)'); ylabel('Reaction time (s)');
changeAxesFontSize(gca,20,20);
export_fig('dots_RT','-eps');


figure(108); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
if useGauss
    beta = [amplConf muConf sigmaConf baselineConf];
    h = plot(xVals, flippedGauss(beta,xVals),'r-','LineWidth',2); hold on;
    errorbar(cohs, pHigh, pHighSE, 'o', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
else
    errorbar(cohs, pHigh, pHighSE, 'o-', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
end
set(gca,'xtick',-0.5:0.25:0.50,'tickdir','out','box','off');
ylim([0.5 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (% coh)'); ylabel('Proportion high bet');
changeAxesFontSize(gca,20,20);
export_fig('dots_conf','-eps');

figure(109); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
errorbar(cohs, pHighCorr, pHighSEcorr, 'o-', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'LineWidth', 2); hold on;
errorbar(cohs, pHighErr, pHighSEerr, 'o--', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
set(gca,'xtick',-0.5:0.25:0.50,'tickdir','out','box','off');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('Motion strength (% coh)'); ylabel('Proportion high bet');
h = legend('Correct','Error','Location','northeast');
set(h,'FontSize',14); legend('boxoff');
set(h,'Position',[0.6595    0.8328    0.2444    0.1203]);
changeAxesFontSize(gca,20,20);
export_fig('dots_conf_split','-eps');







