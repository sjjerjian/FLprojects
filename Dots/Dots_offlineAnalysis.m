% offline analysis wrapper for PLDAPS data, Dots protocols
% CF spawned it from dots3DMP_offlineAnalysis, begun 10-3-18

% skip this cell to analyze a data struct already in the workspace!

clear all; close all

% these will specify which (previously saved) mat file to load
subject = 'hanzo';

% all Hanzo good data
% dateRange = 20190401:20200315; % good (PDW, not RT)
% OR
% dateRange = 20200101:20200731; % large rightward bias on low bets, where
                                 % does it come from?
% dateRange = 20200101:20200315; % same. gotta go month by month

% Hanzo for R01
% dateRange = 20200728:20200930; % great PDW; misnamed (not really up to 9/30)
dateRange = 20200820:20200922; % really nice RT
% dateRange = 20200901:20200922; % good for both, with TrNo<600 or even 800

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
% should really use a two-sample Z test for proportions...



%% plot

Dots_plot

% for nicer looking graphs:
% Dots_plot_forTalk


%% fit simple DDM [code merge in progress]

% fitDDM_simple(data.scoh,data.choice-1,data.RT); % choice-1 because coded as 1=left, 2=right




