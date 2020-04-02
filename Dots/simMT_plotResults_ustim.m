

% calculate prob rightward, high bet
pRight = nan(length(cohs),2);
pHigh = nan(length(cohs),2);
RTmean = nan(length(cohs),2);
RTse = nan(length(cohs),2);
nc = nan(length(cohs),2);
for c = 1:length(cohs)
    for s = 0:1
        I = coh==cohs(c) & ustim==s;
        pRight(c,s+1) = sum(I & choice==1) /  sum(I);
        pHigh(c,s+1) = sum(I & pdw==1) /  sum(I);
        RTmean(c,s+1) = mean(RT(I));
        RTse(c,s+1) = std(RT(I)) / sqrt(sum(I));
        nc(c,s+1) = sum(I);
    end
end
pRightSE = sqrt(pRight.*(1-pRight)./nc);
pHighSE = sqrt(pHigh.*(1-pHigh)./nc);

figure;  clf; set(gcf,'Position',[900 320 480 920]);

% choice
subplot(3,1,1);
errorbar(cohs,pRight(:,1),pRightSE(:,1),'bo'); hold on;
errorbar(cohs,pRight(:,2),pRightSE(:,2),'ro');
% logistic fit
xVals = min(coh):0.01:max(coh);
[beta, ~, ~] = glmfit(coh(ustim==0), choice(ustim==0), 'binomial');
yVals(:,1) = glmval(beta,xVals,'logit');
plot(xVals,yVals,'b');
[beta, ~, ~] = glmfit(coh(ustim==1), choice(ustim==1), 'binomial');
yVals(:,1) = glmval(beta,xVals,'logit');
plot(xVals,yVals,'r');
ylabel('P(right choice)');
legend('no ustim','ustim');

% RT
subplot(3,1,2);
errorbar(cohs,RTmean(:,1),RTse(:,1),'bs-'); hold on;
errorbar(cohs,RTmean(:,2),RTse(:,2),'rs-');
ylabel('RT (ms)'); % again note this isn't realistic for variable dur

% PDW
subplot(3,1,3);
errorbar(cohs,pHigh(:,1),pHighSE(:,1),'bo-'); hold on;
errorbar(cohs,pHigh(:,2),pHighSE(:,2),'ro-');
ylim([0 1]); xlabel('coh'); ylabel('P(high bet)');