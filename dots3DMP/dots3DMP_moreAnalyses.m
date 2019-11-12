%% does overall conf depend on delta? (nope)

exportfigs = 0;

for c = 1:length(cohs)
    I = data.modality==3 & data.coherence==cohs(c) & data.delta==0;
    confDzero(c) = mean(data.conf(I));
    confDzeroSE(c) = std(data.conf(I))/sqrt(sum(I));

    J = data.modality==3 & data.coherence==cohs(c) & data.delta~=0;
    confDnonzero(c) = mean(data.conf(J));
    confDnonzeroSE(c) = std(data.conf(J))/sqrt(sum(J));

    [~,pval(c)] = ttest2(data.conf(I),data.conf(J))
end

barx = 1:length(cohs);
bary = [confDzero ; confDnonzero]';
errlow = [confDzeroSE ; confDnonzeroSE]';
errhigh = errlow;
figure(809); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
bar(barx,bary); hold on;
er = errorbar([barx-0.14;barx+0.14]',bary,errlow,errhigh,'k.');    
ylim([0 1]);
set(gca,'XTickLabel',[10 50 90],'Ytick',0:0.2:1);
xlabel('Visual coherence (%)'); ylabel('Avg confidence rating');
legend('no conflict', 'conflict'); legend('boxoff');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
if exportfigs; export_fig('confDeltaZeroNonzero','-eps'); end


%% calculate predicted/empirical weights

for c = 1:length(cohs)      % m c d
    wvesPred(c) = (1/sigmaPMF(1,1,D)^2) / ((1/sigmaPMF(1,1,D)^2) + (1/sigmaPMF(2,c,D)^2));
                     % m c d
    actual(1) = (muPMF(3,c,1)-muPMF(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (muPMF(3,c,3)-muPMF(3,c,2)+(deltas(3)/2)) / deltas(3);    
    wvesEmp(c) = mean(actual);

    actual(1) = (muConf(3,c,1)-muConf(3,c,2)+(deltas(1)/2)) / deltas(1);
    actual(2) = (muConf(3,c,3)-muConf(3,c,2)+(deltas(3)/2)) / deltas(3);        
    wvesConfBased(c) = mean(actual);
end

nboots = 250;
tic
dots3DMP_cgauss_bootstrap
toc

% looky there!
figure(810); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
er = errorbar(cohs,wvesPred,std(wvesPredboot),std(wvesPredboot),'k-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','k'); hold on;
set(er,'Color','k');
ylim([0 1]);
set(gca,'XTick',cohs,'XTickLabel',[10 50 90],'Ytick',0:0.2:1);
xlabel('Visual coherence (%)'); ylabel('Vestibular weight');
changeAxesFontSize(gca,20,20); set(gca,'box','off');
l = legend('Pred from single-cues'); legend('boxoff');
set(l,'Position',[0.3312    0.7422    0.6069    0.2172]);
if exportfigs; export_fig('weights1','-eps'); end

er = errorbar(cohs,wvesConfBased,std(wvesConfBasedboot),std(wvesConfBasedboot),'c-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','c');    
set(er,'Color','c')
l = legend('Pred from single-cues', 'Pred from conf');
set(l,'Position',[0.3312    0.7422    0.6069    0.2172]);
if exportfigs; export_fig('weights2','-eps'); end

er = errorbar(cohs,wvesEmp,std(wvesEmpboot),std(wvesEmpboot),'m-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','m');    
set(er,'Color','m')
l = legend('Pred from single-cues', 'Pred from conf','Empirical');
set(l,'Position',[0.3312    0.7422    0.6069    0.2172]);
if exportfigs; export_fig('weights3','-eps'); end




%% compare with conf shifts:

% first just raw biases
n = 1;
for c=1:length(cohs)
    for d = 1:length(deltas)
        choiceBias(n) = muPMF(3,c,d);
        choiceBiasSE(n) = muPMFse(3,c,d);
        confBias(n) = muConf(3,c,d);
        confBiasSE(n) = muConfse(3,c,d);
        n = n+1;
    end
end
figure(20); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
[hsym,hxe,hye] = errorbar2(choiceBias, confBias, choiceBias, confBias, 'ko', 2);
set(hsym,'MarkerSize',10,'MarkerFaceColor','w');
hold on; plot([-2 2],[-2 2],'k--','LineWidth',2); axis square;
xlim([-2 2]); ylim([-2 2]);
set(gca,'Xtick',-2:1:2,'Ytick',-2:1:2);
xlabel('Choice shift'); ylabel('Confidence shift');
changeAxesFontSize(gca,20,20); set(gca,'box','off')
if exportfigs; export_fig('biasVsConf','-eps'); end


%% relationship bw conf and weights

% for each coh, bin trials by conf (try median split first), then calc
% weights separately for low and high conf trials

% this logic doesn't work, because confidence in the combined condition per
% se would not map onto weights even under the hypothesis

% BUT we can do it peacewise:

% high coh, left side of the graph, delta +3,
I = data.modality==3 & ismember(data.coherence,cohs(2:end)) & data.heading<=hdgs(3) & data.delta==deltas(end);
% when conf is lower, Pright should be higher (anticorrelated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

clr = {[255 180 0]./255, [175 238 238]./255, [255 140 255]./255};
figure(21); set(gcf,'Color',[1 1 1],'Position',[50 20 250 200],'PaperPositionMode','auto'); clf;
h = plot([1 2], [PrightLow PrightHigh], 'ok-', 'MarkerSize',10, 'MarkerFaceColor','k');
set(h,'Color', clr{3}, 'MarkerFaceColor',clr{3});
xlim([0.5 2.5]); ylim([0 1]);
set(gca,'Xtick',[1 2],'XTickLabel',{'Low Conf','High Conf'},'Ytick',0:0.2:1);
ylabel('Proportion Rightward');
changeAxesFontSize(gca,16,16); set(gca,'box','off');
hold on;
if exportfigs; export_fig('medianSplit1','-eps'); end 

%% low coh, left side of the graph, delta -3, 
I = data.modality==3 & data.coherence==cohs(1) & data.heading<=hdgs(3) & data.delta==deltas(1);
% when conf lower, Pright should be higher (anticorrelated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'o-', 'MarkerSize', 10);
set(h,'Color', clr{3}, 'MarkerFaceColor', clr{3});
if exportfigs; export_fig('medianSplit2','-eps'); end


%% high coh, right side of the graph, delta -3, 
I = data.modality==3 & ismember(data.coherence,cohs(2:end)) & data.heading>=hdgs(end-2) & data.delta==deltas(1);
% when conf is lower, Pright should be lower (correlated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'ob-', 'MarkerSize',10, 'MarkerFaceColor','b');
set(h,'Color', clr{1}, 'MarkerFaceColor',clr{1});
if exportfigs; export_fig('medianSplit3','-eps'); end


%% low coh, right side of the graph, delta +3, 
I = data.modality==3 & data.coherence==cohs(1) & data.heading>=hdgs(end-2) & data.delta==deltas(end);
% when conf lower, Pright should be lower (correlated)
confLow = data.conf(I) < median(data.conf(I));
confHigh = data.conf(I) >= median(data.conf(I));
PrightLow = sum(confLow & data.choice(I)==2) / sum(confLow);
PrightHigh = sum(confHigh & data.choice(I)==2) / sum(confHigh);

h = plot([1 2], [PrightLow PrightHigh], 'ob-', 'MarkerSize',10, 'MarkerFaceColor','b');
set(h,'Color', clr{1}, 'MarkerFaceColor',clr{1});
if exportfigs; export_fig('medianSplit4','-eps'); end





%% try sliding window across trials (or session by session!)






%% for RT, try kiani 2014 analysis:

% plot conf as a function of RT quantile, separately for abs(hdg)

clear X Y
uhdg = unique(abs(data.heading));
for h = 1:length(uhdg)
    I = abs(data.heading)==uhdg(h);
    theseRT = data.rt(I);
    theseConf = data.conf(I);
    rtQ = [0 quantile(theseRT,3) inf]; % try quartiles
    for q = 1:length(rtQ)-1
        J = theseRT>=rtQ(q) & theseRT<rtQ(q+1);
        X(h,q) = mean(theseRT(J));
        Y(h,q) = mean(theseConf(J));
    end
end

Ltxt = {num2str(uhdg(1)), num2str(uhdg(2)), num2str(uhdg(3)), num2str(uhdg(4)), num2str(uhdg(5))};
figure(16); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
clr = {'bo-','rs-','c^-','gd-','kv-'};
clear g L
for h = 1:length(uhdg)       % vvvv  shifted for clarity!!!
    g(h) = plot(X(h,:),Y(h,:)+(0.03*double(h==5)),clr{h},'LineWidth', 2); hold on;
    set(g(h),'MarkerSize',10,'MarkerFaceColor',clr{h}(1));
    xlim([1 3.5]); ylim([0 1]);
    % set(gca,'Xtick',-2:1:2,'Ytick',-2:1:2);
    xlabel('Reaction time (s)'); ylabel('Confidence');
    changeAxesFontSize(gca,16,16); set(gca,'box','off');
    L{h} = Ltxt{h};
    l = legend(g,L,'Location','South','Orientation','Horizontal');
    legend('boxoff')
    export_fig(['kianiFig2-' num2str(h)],'-eps');
end






%% fit DDM

options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';

% params: 

% % % 
% % % 
% % % % Drugowitsch model has 12 params, plus 8 for biases and lapse rates
% % % % we'll skip the latter, and we can drop the 3 Tnd terms until we are fitting RT
% % % % so we have 9 params:
% % % 
% % %     %    aVis gammaVis bVis thetaVis kVes thetaVes gammaCom bCom thetaCom
% % % fixed = [0    0        0    0        0    0        0        0    0       ];
% % % 
% % % % per Drugowitsch, variance of momentary evidence scales with coherence as:
% % % % var(c) ~ 1 + bVis*coh^gammaVis;
% % % 
% % % % similarly,
% % % % kVis ~ aVis*cVis^gammaVis
% % % 
% % % % the bound gets a free parameter (theta) for each modality, but this is
% % % % only because the variance is what really changes across conditions, and
% % % % this is absorbed into the definition of (normalized) bounds. I'm still
% % % % not sure about this...
% % %  
% % % % initial guess (or hand-tuned params)
% % % aVis = 0.5; % sensitivity parameter, multiplies coh
% % % gammaVis = 1; % determines scaling of sensitivity (and variance) by coh
% % % bVis = 30; % bound height
% % % theta = 1.6; % criterion (in log odds correct) for betting high
% % % alpha = 0.1; % base rate of low-bet choices
% % % 
% % % guess = [aVis gammaVis bVis thetaVis kVes thetaVes gammaCom bCom thetaCom];
% % % 
% % % 
% % % 
% % % 


% put this on hold for now, instead use a much simpler version,
% the minimum parameterization sufficient for CCC_sim


data.strength = sind(data.heading);
data.dur = ones(size(data.strength))*2000;

    %    kves kvis B 
fixed = [0    0    0];

% one small diff: in sim, kvis is just coh, here it will multiply coh

% initial guess (or hand-tuned params)
kves = 0.4; % sensitivity parameter for ves
kvis = 1; % multiplies coh to get sensitivity parameter for vis
B = 70; % bound height

guess = [kves kvis B];


% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=1;
% ************************************

options.feedback = false;
options.plot = false;

[X, LL_final, data, fit] = dots3DMP_fitDDM(data,options,guess,fixed);





%% calculate 'score' based on accuracy and metacog accuracy

% what we don't want is simply to reward confidence scaling with heading
% angle, or coherence
% rather, within a coh+hdg, how well does confidence predict accuracy?

corr = double((data.heading>0 & data.choice==2) | (data.heading<0 & data.choice==1));
corr(data.heading==0) = NaN;
q=1;
for m = 1:length(mods)
    for c = 1:length(cohs)
        for h = 1:length(hdgs)
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h); % all trials irrespective of delta
            confCorr = mean(data.conf(J & corr==1));
            confErr = mean(data.conf(J & corr==0));
%             diff(m,c,h) = confCorr - confErr;
            diff(q) = confCorr - confErr;

            cc = corrcoef(corr(J),data.conf(J));
%             phi(m,c,h) = cc(1,2);
            phi(q) = cc(1,2);
            
            q = q+1;
        end
    end
end



meanphi = nanmean(phi);
bestphi = 0.25; % hypothetical best performance
worstphi = 0; % not sure about this

totalPcorr = nansum(corr)/sum(data.heading~=0);
bestPcorr = 0.9; % hypothetical best performance
worstPcorr = 0.5;

scorePcorr = (min([totalPcorr bestPcorr]) - worstPcorr) / (bestPcorr-worstPcorr);
scorePhi = (min([meanphi bestphi]) - worstphi) / (bestphi-worstphi);

score =  scorePcorr*0.5 + scorePhi*0.5
payout = score*8 + 12




