

% calculate prob rightward, high bet
pRight = nan(length(cohs),1);
pHigh = nan(length(cohs),1);
RTmean = nan(length(cohs),1);
RTse = nan(length(cohs),1);
nc = nan(length(cohs),1);
for c = 1:length(cohs)
    I = coh==cohs(c);
    pRight(c,1) = sum(I & choice==1) /  sum(I);
    pHigh(c,1) = sum(I & pdw==1) /  sum(I);
    RTmean(c,1) = mean(RT(I));
    RTse(c,1) = std(RT(I)) / sqrt(sum(I));
    nc(c,1) = sum(I);
%     meanLogOddsCorr(c,1) = mean(logOddsCorr(I)); % temp
end
pRightSE = sqrt(pRight.*(1-pRight)./nc);
pHighSE = sqrt(pHigh.*(1-pHigh)./nc);


% % temp: diagnose mismatch between logOddsCorr and actual pcorr
% pCorrActual = pRight;
% pCorrActual(cohs<0) = 1-pCorrActual(cohs<0)
% pCorrFromLogOdds = logistic(meanLogOddsCorr)


figure;  clf; set(gcf,'Position',[900 320 480 920]);
% choice
subplot(3,1,1); errorbar(cohs,pRight,pRightSE,'o'); hold on;
% logistic fit
[beta, ~, ~] = glmfit(coh, choice, 'binomial');
xVals = min(coh):0.01:max(coh);
yVals = glmval(beta,xVals,'logit');
plot(xVals,yVals);
ylabel('P(right choice)');
% RT
subplot(3,1,2);
errorbar(cohs,RTmean,RTse,'s-');
ylabel('RT (ms)'); % again note this isn't realistic for variable dur
% PDW
subplot(3,1,3);
errorbar(cohs,pHigh,pHighSE,'o-');
ylim([0 1]); xlabel('coh'); ylabel('P(high bet)');