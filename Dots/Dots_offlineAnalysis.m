% offline analysis wrapper for PLDAPS data, Dots protocols
% CF spawned it from dots3DMP_offlineAnalysis, begun 10-3-18

% skip this cell to analyze a data struct already in the workspace!

clear all; close all

% these will specify which (previously saved) mat file to load
subject = 'hanzo';

% all Hanzo good data
% dateRange = 20190401:20201231; % all 'good' data

% Hanzo for R01
% dateRange = 20200728:20200930; % great PDW; misnamed (not really up to 9/30)
% dateRange = 20200820:20200922; % nice RT, but PDW kinda flat
% dateRange = 20200901:20200922; % good for both, with TrNo<600 or even 800

% some indiv months
dateRange = 20190501:20190531; % var dur example
% dateRange = 20201101:20201130; % RT example

folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
load([folder file], 'data');


%% more cleanup

% clean up RT variable for mixed RT/non-RT datasets
if isfield(data,'RT')
    data.RT(data.RT>1.5)=nan; % anything over 1.5s is almost certainly invalid,
                              % and this also gets rid of the infs
    data.RT(data.RT<0.01)=nan; % zeroes too
else
    data.RT = nan(size(data.choice));
end

% not doing anyhting with these for now:
try
    data = rmfield(data,{'forgivePDW','trainingSeqPDW','confidenceTime'});
catch
end

% remove invalid trials (fixation breaks, etc)
% removethese = isnan(data.choice) | isnan(data.PDW); % may need to add conditions here
removethese = isnan(data.choice); % not all sessions have PDW! nans will be removed below
% removethese = isnan(data.choice) | ~isnan(data.PDW) % and for RT only, may want to exclude PDW trials
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% define signed coherence; positive = rightward (0 deg)
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);

% random fixes
data.scoh(abs(data.scoh)<0.001)=0;


% % quick look at blocks, for when some need to be excluded
% blocks = unique(data.filename);
% nTrialsByBlock = nan(length(blocks),1);
% data.trialNum = zeros(size(data.filename));
% for u = 1:length(blocks)
%     nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
%     % also number trials within a block, to look at early vs late trials
%     data.trialNum(strcmp(blocks(u),data.filename)) = 1:sum(strcmp(blocks(u),data.filename))';
% end

% we can be pretty sure blocks with <N trials (say, 10) are to be discarded
% removethese = ismember(data.filename,blocks(nTrialsByBlock<10));


% % remove early/late trials
% removethese = removethese | data.trialNum>800;
% fnames = fieldnames(data);
% for F = 1:length(fnames)
%     eval(['data.' fnames{F} '(removethese) = [];']);
% end



%% parse data

Dots_parse


%% some basic stats

% t test on slopes (is sensitivity greater on high-bet trials?):
mu = abs(B2(2)-B3(2));
se = sqrt(stats2.se(2)^2 + stats3.se(2)^2);
t = mu/se;
df = sum(~isnan(data.PDW))-length(B2); 
pval_ttest_slopes = 2*(1-tcdf(t,df)) % two-tailed


% t test on pHigh, corr vs err (is confidence higher on correct trials?)
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
% [should really use a two-sample Z test for proportions?]



%% plot

Dots_plot


% % for nicer looking graphs:
% Dots_plot_forTalk


%% if var dur, check conf/accuracy vs. dur
if sum(isnan(data.RT))>0.8*length(data.RT) % arbitrary minimum proportion of non-RT trials; eventually make it a flag
    Dots_TDA;
end



%% fit simple DDM

% for RT data, can start here (ignores PDW, but a decent tutorial, based on
% Shadlen et al. 2006 and Palmer et al. 2005)
if sum(isnan(data.RT))<0.8*length(data.RT) % arbitrary minimum proportion of RT trials; eventually make it a flag
    b = fitDDM_simple(data.scoh,data.choice-1,data.RT); % (choice-1 because coded as 1=left, 2=right)
end


%% for var-dur + PDW, switch to Kiani 09 method

if sum(isnan(data.RT))>0.8*length(data.RT) % arbitrary minimum proportion of non-RT trials; eventually make it a flag

    % rename some vars:
    D.strength = data.scoh;
    D.dur = round(data.duration*1000); % must be integers, in ms
    D.choice = data.choice-1; 
    D.conf = data.PDW;

    options.fitMethod = 'fms'; % fminsearch
    % options.fitMethod = 'global';
    % options.fitMethod = 'multi';
    % options.fitMethod = 'pattern';

    % params: 
    %        1 2 3
    %       [k B theta alpha]
    fixed = [0 0 0 0]; % can fix some params and fit the others

    % % initial guess (or hand-tuned params)
    k = 0.6; % sensitivity parameter
    B = 15; % bound height
    theta = 1.2; % criterion (in log odds correct) for betting high
    alpha = 0.2; % base rate of low-bet choices

    guess = [k B theta alpha];

    % ************************************
    % set all fixed to 1 for hand-tuning:
    fixed(:)=1;
    % ************************************

    options.feedback = 1; % 1 = text output to cmd window, 2 = that and plot LL across runs
    options.plot = 0; % plot the marginal PDFs, logOddsCorr map, and high/low bet regions (only makes sense for fixed(:)=1)

    % fit it!
    [X, LL_final, D, fit] = fitDDM_wConfidence_simple(D,options,guess,fixed);
        
    % plot the fits and compare to data
    cohs_fit = unique(fit.strength);
    pRight_model = nan(length(cohs_fit),1);
    pRightHigh_model = nan(length(cohs_fit),1);
    pRightLow_model = nan(length(cohs_fit),1);
    pHigh_model = nan(length(cohs_fit),1);
    for c = 1:length(cohs_fit)
        I = fit.strength==cohs_fit(c);
        pRight_model(c) = mean(fit.expectedPright(I));
        pRightHigh_model(c) = mean(fit.expectedPrightHigh(I));
        pRightLow_model(c) = mean(fit.expectedPrightLow(I));
        pHigh_model(c) = mean(fit.expectedPhigh(I));
    end
    
    figure; set(gcf,'Position',[86 925 1070 420]);
    subplot(1,2,1);
%     errorbar(cohs,pRight,pRightSE,'bo');
%     hold on; plot(cohs_fit,pRight_model,'c--');
    errorbar(cohs, pRightHigh, pRightSEhigh, 'bo', 'MarkerFaceColor', 'b'); hold on;
    errorbar(cohs, pRightLow, pRightSElow, 'bo', 'MarkerFaceColor', 'w');
    plot(cohs_fit,pRightHigh_model,'c-');
    plot(cohs_fit,pRightLow_model,'c--');
    xlabel('motion strength (%coh)');
    ylabel('proportion rightward choices');

    subplot(1,2,2); errorbar(cohs,pHigh,pHighSE,'ro');
    hold on; plot(cohs_fit,pHigh_model,'m--'); ylim([0 1]);
    xlabel('motion strength (%coh)');
    ylabel('proportion high bets');
    
end


%% for fitting RT+PDW, need to switch to 2D (anticorrelated race) model, e.g. Kiani et al. 2014
if sum(isnan(data.RT))<0.8*length(data.RT) % arbitrary minimum proportion of RT trials; eventually make it a flag

    % [to be merged from dots3DMP track]

end







