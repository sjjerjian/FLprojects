% generic offline analysis wrapper for PLDAPS data, dots3DMP protocol
% CF started it 10-3-18

% this is for the human version with vertical color-gradient targets
% (single-interval, 2AFC heading discrim task, where choice, RT and
% confidence are indicated in a single saccade)

    % struct data has the following fields:
% heading (in degrees relative to straight ahead, positive is rightward)
% choice (1 is left, 2 is right)
% rt (in seconds)
% conf (saccadic endpoint as a fraction of the bar length
    % values >1 or <0 are possible because some overshoot of the bar is 
    % tolerated; could rescale to 0..1 or 0.5..1, with or without
    % rectifying at 0/1 first.


clear all
load devin5days.mat
removethese = isnan(data.choice); % fixation breaks, mostly
data.heading(removethese) = [];
data.choice(removethese) = [];
data.rt(removethese) = [];
data.conf(removethese) = [];

ucoh = unique(data.heading);
pRight = nan(length(ucoh),1);
pRight_se = nan(length(ucoh),1);
RTmean = nan(length(ucoh),1);
RTse = nan(length(ucoh),1);
confMean = nan(length(ucoh),1);
confSE = nan(length(ucoh),1);
N = nan(length(ucoh),1);
for c = 1:length(ucoh)
    I = data.heading==ucoh(c);
    pRight(c) = sum(data.choice(I)==2) / sum(I); % 2 is rightward
    pRight_se(c) = sqrt(pRight(c)*(1-pRight(c)) / sum(I)); % formula for standard error of a proportion
    RTmean(c) = mean(data.rt(I));
    RTse(c) = std(data.rt(I))/sqrt(sum(I));
    confMean(c) = mean(data.conf(I));
    confSE(c) = std(data.conf(I))/sqrt(sum(I));
    N(c) = sum(I);
end

% plot choice and fit logistic regression
% B(1) is the bias, B(2) is the slope, or sensitivity (related to accuracy)
X = data.heading;
y = data.choice==2;
[B, ~, ~] = glmfit(X, y, 'binomial');
xVals = ucoh(1):0.01:ucoh(end);
yVals = glmval(B,xVals,'logit');
figure; errorbar(ucoh,pRight,pRight_se,'bo'); hold on;
plot(xVals,yVals,'-');
xlabel('heading angle (deg)');
ylabel('proportion rightward choices');

% plot RT
figure; errorbar(ucoh,RTmean,RTse,'ro-');
xlabel('signed motion strength (%coh)');
ylabel('reaction time (s)');

% plot conf
figure; errorbar(ucoh,confMean,confSE,'go-');
xlabel('signed motion strength (%coh)');
ylabel('saccadic endpoint (conf) (%)');



