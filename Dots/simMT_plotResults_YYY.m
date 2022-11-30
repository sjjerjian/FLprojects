
% calculate prob rightward, high bet
pRight = nan(length(cohs),1);
pHigh = nan(length(cohs),1);
RTmean = nan(length(cohs),1);
RTse = nan(length(cohs),1);
nc = nan(length(cohs),1);
for c = 1:length(cohs)
    I = coh(1:tr-1)==cohs(c);   % I = logical index to trials of coh 'c'
    pRight(c,1) = sum(I & choice(1:tr-1)==1) /  sum(I);  % count number of trials for choice Right and coh 'c' / total number of 'c' trials
    pHigh(c,1) = sum(I & pdw(1:tr-1)==1) /  sum(I); % Same but for Lo choice trials
    RTmean(c,1) = mean(RT(I)); % mean RT of 'c' trials
    RTse(c,1) = std(RT(I)) / sqrt(sum(I)); % Std error of RT on 'c' trials
    nc(c,1) = sum(I); % store number of trials of coh 'c'
%     meanLogOddsCorr(c,1) = mean(logOddsCorr(I)); % temp
end
pRightSE = sqrt(pRight.*(1-pRight)./nc); % Std error of R choice prob at each coh level
pHighSE = sqrt(pHigh.*(1-pHigh)./nc);

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
errorbar(cohs,RTmean,RTse,'s-'); %RTmean is arranged such that each row is the mean RT of the associated column of cohs
ylabel('RT (ms)'); % again note this isn't realistic for variable dur
% PDW
subplot(3,1,3);
errorbar(cohs,pHigh,pHighSE,'o-');
ylim([0 1]); xlabel('coh'); ylabel('P(high bet)');