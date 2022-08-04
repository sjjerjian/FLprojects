                                  
% load data

clear all

load human_bandit_round1_signed.mat
left_trials = (data.direction == 180);

%load banditDataParseForDist/banditData/round1/humanS5round1.mat
%left_trials = (data.accuracy== 1 & data.guess ==1) | (data.accuracy== 0 & data.guess ==2);
%data.direction = zeros(size(data.accuracy,1), size(data.accuracy,2));
%data.direction(left_trials)=180;

validTrials = ~(data.rt==0);
data.coherence(left_trials) = -data.coherence(left_trials);

% rename some vars:
D.strength = data.coherence(validTrials);
D.dur = round(data.rt(validTrials)*1000);
D.dotchoice = data.guess(validTrials); 
D.dotchoice = D.dotchoice-1;
D.correct = data.accuracy(validTrials);
D.conf = data.conf(validTrials);


options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';


% params: 
%        1 2 3
%       [k B theta alpha]
fixed = [0 0 0 0];


% initial guess (or hand-tuned params)
k = 0.25; % sensitivity parameter
B = 25; % bound height
theta = 1.6; % criterion (in log odds correct) for betting high
alpha = 0.1; % base rate of low-bet choices


guess = [k B theta alpha];

% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=0;
% ************************************

options.feedback = false;
options.plot = false;

[X, LL_final, data, fit] = BanditDots_fitDDM(D,options,guess,fixed);


% plot data and compare to fits

ucoh = unique(data.strength);
pRight = nan(length(ucoh),1);
pRight_se = nan(length(ucoh),1);
pHigh = nan(length(ucoh),1); % high bet (conf)
pHigh_se = nan(length(ucoh),1);
N = nan(length(ucoh),1);
    % SEs of the expected proportions from the model are probably not
    % meaningful. What matters are the SEs of the parameters themselves
for c = 1:length(ucoh)
    I = data.strength==ucoh(c);
    pRight(c) = sum(data.dotchoice(I)==1) / sum(I); % 1 is rightward
    pRight_se(c) = sqrt(pRight(c)*(1-pRight(c)) / sum(I)); % formula for standard error of a proportion
    %pHigh(c) = sum(data.conf(I)>=0.5) / sum(I); % 1 is high bet
    pHigh(c) = mean(data.conf(I));
    pHigh_se(c) = sqrt(pHigh(c)*(1-pHigh(c)) / sum(I)); 
    N(c) = sum(I);
end

% now the model
ucoh_fit = unique(fit.strength);
pRight_model = nan(length(ucoh_fit),1);
pHigh_model = nan(length(ucoh_fit),1); % high bet (conf)
for c = 1:length(ucoh_fit)
    I = fit.strength==ucoh_fit(c);
    pRight_model(c) = mean(fit.expectedPright(I));
    pHigh_model(c) = mean(fit.expectedPhigh(I));
        % should separate this out by duration for a more fine-grained look.
        % Then there's also confidence conditioned on correct/error, and
        % accuracy conditioned on high/low bet...
end

figure; set(gcf,'Position',[86 925 1070 420]);
subplot(1,2,1); errorbar(ucoh,pRight,pRight_se,'bo-', 'LineWidth', 1.5);
hold on; plot(ucoh_fit,pRight_model, '--', 'LineWidth', 1.5);
xlabel('motion strength (%coh)');
ylabel('proportion rightward choices');

subplot(1,2,2); errorbar(ucoh,pHigh,pHigh_se,'ro-', 'LineWidth', 1.5);
hold on; plot(ucoh_fit,pHigh_model, '--', 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]); ylim([0 1]);
xlabel('motion strength (%coh)');
ylabel('proportion high bets');

%% Plot the log odds correct and the probability density map

% dots task params
cohs = unique(data.strength);
dT = 1; % ms
maxdur = max(max(data.dur)); % max viewing duration (RT task)
timeAxis = 0:dT:maxdur;
sigma = 1; % SD of momentary evidence
theta = 1; % threshold for Psure choice in units of log odds correct (Kiani & Shadlen 2009)

% Kiani map [need to get 2014 version!]
[logOddsMapR, logOddsMapL, logOddsCorrMap, PMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth_hy(k,B,sigma,theta,cohs,timeAxis,0);


% define the grid resolution for time and decision variable
dt = 0.1; % 0.1 ms seems to work best for FP4 even though
% we don't store the vars at this resolution
dx = min(0.1, B/100);
% define the time axis 
max_dur = max(data.dur);
t = dt:dt:max_dur;
% define xmesh and the ghost cells for bound crossing probability
b_margin = repmat(4*sigma, [1 2]);
xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';

n=30; % number of color levels (smoothness)
figure; x = repmat(1:size(logOddsCorrMap, 2),size(logOddsCorrMap, 1),1)';
y = repmat(vAxis',1,size(logOddsCorrMap, 2))';

subplot(2,1,1);         

[~,h] = contourf(x,y,logOddsCorrMap',n); title('Log odds correct'); colorbar;
if n>20; set(h,'LineColor','none'); end
xlabel('t (ms)'); ylabel('x');
        
subplot(2,1,2); 

[~,h] = contourf(x,y,PMap',n); title('Probability density'); colorbar;
if n>20; set(h,'LineColor','none'); end
xlabel('t (ms)'); ylabel('x');
caxis([0 0.004])
        
