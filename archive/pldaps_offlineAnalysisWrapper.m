% generic offline analysis wrapper for PLDAPS data
% CF started it 10-3-18

clearvars -except file

% decide which files to load
% folder = '/home/fetschlab/data/';
subject = 'hanzo';
localDir = ['/Users/chris/Documents/MATLAB/PLDAPS_data/' subject '/'];
dateRange = 20181119:20181121;
protocol = 'DotsBasic';

% set up data structure
data = struct;
data.dir = [];
data.coh = [];
% data.heading = [];
data.choice = [];
data.rt = [];
data.conf = [];
    % add variables here

% search folder for matching files and extract the desired variables from PDS
allFiles = dir(localDir);
for d = 1:length(dateRange)
    for f = 3:length(allFiles) % skip 1+2, they are are "." and ".."
        if contains(allFiles(f).name, subject) ... % check if the target file is in localDir
           && contains(allFiles(f).name, num2str(dateRange(d))) ...
           && contains(allFiles(f).name, protocol)
       
            disp(['loading ' allFiles(f).name]);
            load([localDir allFiles(f).name],'-mat'); % load it. this is the time limiting step;
                                                      % will eventually change how data are saved to make this faster
                                                      
            for t = 1:length(PDS.data) % loop over trials for this file
                if isfield(PDS.data{t}.behavior,'choice') % and save out the data, excluding trials
                                                          % with missing data for whatever reason
                    data.dir(end+1,1) = PDS.conditions{t}.stimulus.direction;
                    data.coh(end+1,1) = PDS.conditions{t}.stimulus.coherence;
%                     data.heading(end+1,1) = PDS.conditions{t}.stimulus.heading;
                    data.choice(end+1,1) = PDS.data{t}.behavior.choice;
                    try data.rt(end+1,1) = PDS.data{t}.behavior.RT; catch isRTtask = 0; end
                    try data.conf(end+1,1) = PDS.data{t}.behavior.saccEndPoint;
                end
            end
            clear PDS
        end
    end
end

% save data so you don't have to repeat the time consuming step
file = [subject '_' num2str(dateRange(1)-20180000) '_' num2str(dateRange(end)-20180000) '.mat'];
localDir = '/Users/chris/Documents/MATLAB/';
save([localDir file], 'data');


%% plot choice (and RT), and fit weibull (unsigned coh)

clearvars -except file folder
load([localDir file]);

ucoh = unique(data.coh);
nTot = nan(length(ucoh),1);
pCorr = nan(length(ucoh),1);
pCorr_se = nan(length(ucoh),1);
for c = 1:length(ucoh)
    I = data.coh==ucoh(c) & ismember(data.choice,[1 2]);
    % a trial is correct if dir==0 and choice==2 (rightward), OR dir==180 and choice==1 (leftward)
    nCorr = sum( (data.dir(I)==0 & data.choice(I)==2) | (data.dir(I)==180 & data.choice(I)==1) );
    nTot(c) = sum(I);
    pCorr(c) = nCorr / nTot(c);
    pCorr_se(c) = sqrt(pCorr(c)*(1-pCorr(c)) / sum(I)); % formula for standard error of a proportion
end
figure; errorbar(ucoh,pCorr,pCorr_se,'o-');
ylim([0.45 1]);
xlabel('unsigned motion strength (%coh)');
ylabel('proportion rightward choices');

% fit cumulative weibull distribution (aka 'Quick' function)
% alpha is the threshold (coherence level yielding 75% correct), beta is the slope
[alpha, beta, ~, abse, ~] = quickfit([ucoh pCorr nTot]);
xVals = ucoh(1):0.01:ucoh(end);
yVals = 1 - .5 * exp( -(xVals/alpha).^beta ); % the Quick function
figure; errorbar(ucoh,pCorr,pCorr_se,'o'); hold on;
plot(xVals,yVals,'-');
ylim([0.45 1]);
xlabel('unsigned motion strength (%coh)');
ylabel('proportion rightward choices');


%% fit logistic (signed coh)

% convert coherence to signed coherence, using direction (rightward will be positive, leftward negative)
data.scoh = data.coh;
data.scoh(data.dir==180) = data.scoh(data.dir==180) * -1;

uscoh = unique(data.scoh);
pRight = nan(length(uscoh),1);
pRight_se = nan(length(uscoh),1);
for c = 1:length(uscoh)
    I = data.scoh==uscoh(c) & ismember(data.choice,[1 2]);
    pRight(c) = sum(data.choice(I)==2) / sum(I); % 2 is rightward
    pRight_se(c) = sqrt(pRight(c)*(1-pRight(c)) / sum(I)); % formula for standard error of a proportion
end
figure; errorbar(uscoh,pRight,pRight_se,'o-');
xlabel('signed motion strength (%coh)');
ylabel('proportion rightward choices');

% fit logistic regression
% B(1) is the bias, B(2) is the slope, or sensitivity (related to accuracy)
X = data.scoh;
y = data.choice==2;
[B, ~, ~] = glmfit(X, y, 'binomial');
xVals = uscoh(1):0.01:uscoh(end);
yVals = glmval(B,xVals,'logit');
figure; errorbar(uscoh,pRight,pRight_se,'o'); hold on;
plot(xVals,yVals,'-');
xlabel('motion strength (%coh)');
ylabel('proportion rightward choices');


%% function fitDDM_simple(stimStr,choice,RT)

fitDDM_simple(data.scoh,data.choice-1,[]);

