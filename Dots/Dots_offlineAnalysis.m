% offline analysis wrapper for PLDAPS data, Dots protocols
% CF spawned it from dots3DMP_offlineAnalysis, begun 10-3-18

%*****************************
% skip this cell to analyze a data struct already in the workspace!
%*****************************

clear; close all

% these will specify which (previously saved) mat file to load
subject = 'hanzo';

folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';

% all Hanzo good data
% dateRange = 20190401:20201231; % all 'good' data

% Hanzo for R01
% dateRange = 20200728:20200930; % great PDW; misnamed (not really up to 9/30)
% dateRange = 20200820:20200922; % nice RT, but PDW kinda flat
% dateRange = 20200901:20200922; % good for both, with TrNo<600 or even 800

% some indiv months
% dateRange = 20190501:20190531; % var dur example
% dateRange = 20201101:20201130; % RT example

% should be best!
dateRange = 20200801:20201130;

% dateRange = 20210208:20210212; % last week

% maxTrialNum = 700; % set to ~600-800 to omit late trials

if sum(diff(dateRange)>1)==0
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
elseif sum(diff(dateRange)>1)==1
    file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(diff(dateRange)>1)) '+' num2str(dateRange(find(diff(dateRange)>1)+1)) '-' num2str(dateRange(end)) '.mat'];
else
    file = [subject '_' num2str(dateRange(1)) '---' num2str(dateRange(end)) '.mat'];
end

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
removethese = isnan(data.choice);
try removethese = removethese | data.oneTargPDW==1; catch; end 
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% define signed coherence; positive = rightward (0 deg)
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);

% random fixes
data.scoh(abs(data.scoh)<0.001)=0;
if max(data.choice)==2
    data.choice = data.choice-1; % ensure choice is 0 or 1 (sometimes coded as 1..2)
end

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
if exist('maxTrialNum','var')
    removethese = removethese | data.trialNum>maxTrialNum;
    fnames = fieldnames(data);
    for F = 1:length(fnames)
        eval(['data.' fnames{F} '(removethese) = [];']);
    end
end


%% parse data

Dots_parse


% %% some basic stats
% 
% % t test on slopes (is sensitivity greater on high-bet trials?):
% mu = abs(B2(2)-B3(2));
% se = sqrt(stats2.se(2)^2 + stats3.se(2)^2);
% t = mu/se;
% df = sum(~isnan(data.PDW))-length(B2); 
% pval_ttest_slopes = 2*(1-tcdf(t,df)) % two-tailed
% 
% 
% % t test on pHigh, corr vs err (is confidence higher on correct trials?)
% M = ~isnan(data.PDW) & data.correct==1;
% pHighCorr_all = sum(M & data.PDW==1) / sum(M);
% MM = ~isnan(data.PDW) & data.correct==0;
% pHighErr_all = sum(MM & data.PDW==1) / sum(MM);
% 
% pHighSEcorr_all = sqrt( (pHighCorr_all.*(1-pHighCorr_all)) ./ sum(M) );
% pHighSEerr_all = sqrt( (pHighErr_all.*(1-pHighErr_all)) ./ sum(MM) );
% 
% mu = abs(pHighCorr_all-pHighErr_all);
% se = sqrt(pHighSEcorr_all^2 + pHighSEerr_all^2);
% t = mu/se; 
% df = sum(~isnan(data.PDW))-length(B2); 
% pval_ttest_pHighCorrErr = 2*(1-tcdf(t,df)) % two-tailed
% % [should really use a two-sample Z test for proportions?]



%% plot

Dots_plot


% % for nicer looking graphs:
Dots_plot_forTalk


% %% if var dur, check conf/accuracy vs. dur
% if sum(isnan(data.RT))>0.8*length(data.RT) % arbitrary minimum proportion of non-RT trials; eventually make it a flag
%     Dots_TDA; % time-dependent accuracy function
% end



%% fit simplest DDM

% for RT data, can start here (ignores PDW, but a decent tutorial, based on
% Shadlen et al. 2006 and Palmer et al. 2005)
if sum(isnan(data.RT))<0.8*length(data.RT) % arbitrary minimum proportion of RT trials; eventually make it a flag
    [b,~] = fitDDM_simple(data.scoh,data.choice,round(data.RT*1000));
end
% seems to find local minima often -- may need multiple starting points



%% fit DDM with confidence (with or without RT)


%****** first select which model to fit ********
%
% options.errfcn = @errfcn_DDM_1D_wConf; modelID=1; % 1D DDM with threshold on log odds, usually for var dur [Kiani 09 (FP4)]
options.errfcn = @errfcn_DDM_2D_wConf; modelID=2; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (WolpertMOI)]
%
%***********************************************


options.feedback = 1; % 1 = text output to cmd window, 2 = that, plus plot LL across runs
options.plot = 0; % plot the marginal PDFs, logOddsCorr map, and high/low bet regions (only makes sense for fixed(:)=1)

% choose optimization method
options.fitMethod = 'fms'; % fminsearch
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% % options.fitMethod = 'bads'; % not implemented yet, see dots3DMP for work-in-progress

% params: 
switch modelID
    case 1 %errfcn_DDM_1D_wConf
        
        % initial guess (or hand-tuned params)
        k = 0.6; % sensitivity parameter
        B = 15; % bound height
        theta = 1.2; % criterion (in log odds correct) for betting high
        alpha = 0.2; % base rate of low-bet choices

        guess = [k B theta alpha];
        fixed = [0 0 0     0]; % can fix some params and fit the others, or fix all to hand-tune

        data.dur = round(data.duration*1000); % dur must be integer valued (in ms)
        
    case 2 %errfcn_DDM_2D_wConf
        
        k = 7;
        B = 1.5;
        theta = 1;
        alpha = 0.2; % base rate of low-bet choices
        Tnd = 0.3; % non-decision time

        guess = [k B theta alpha Tnd];
        fixed = [0 0 0     0     0];
        
%    case [others yet to be written: extrema, snapshot, SDT models,...]

end

% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=0;
% ************************************

% fit it!
[X, err_final, fit] = Dots_fitDDM_wrapper(guess,fixed,data,options);


%% plot DDM fits + compare to data

% to generate smooth curves, call errfcn again with interpolated
% coh axis (and resampled durs if var dur task)

switch modelID
    case 1 %errfcn_DDM_1D_wConf

        ncohsteps = 33; % should be odd so there's a zero
        g_str = linspace(min(data.scoh),max(data.scoh),ncohsteps);
        ndursamples = 1000;
        rand_dur = randsample(data.dur, ndursamples, 'true');

        fitInterp = struct;
        fitInterp.scoh = repmat(g_str', ndursamples, 1);
        for j = 1:ndursamples
            fitInterp.dur((j-1)*ncohsteps+1 : (j-1)*ncohsteps+ncohsteps) = rand_dur(j);
        end
        fitInterp.dur = fitInterp.dur';
        % just fill these with ones, since error calc doesn't matter, and
        % after this only expectedPright/High is used for the plot
        fitInterp.choice = ones(size(fitInterp.scoh));
        fitInterp.PDW = ones(size(fitInterp.scoh));

        options.plot = 0; options.feedback=0; fixed(:)=1;
        [~,~,fitInterp] = options.errfcn(X,X,fixed,fitInterp,options);
                
        cohs_fit = unique(fitInterp.scoh);
        pRight_model = nan(length(cohs_fit),1);
        pRightHigh_model = nan(length(cohs_fit),1);
        pRightLow_model = nan(length(cohs_fit),1);
        pHigh_model = nan(length(cohs_fit),1);
        for c = 1:length(cohs_fit)
            I = fitInterp.scoh==cohs_fit(c);
            pRight_model(c) = mean(fitInterp.expectedPright(I));
            pRightHigh_model(c) = mean(fitInterp.expectedPrightHigh(I));
            pRightLow_model(c) = mean(fitInterp.expectedPrightLow(I));
            pHigh_model(c) = mean(fitInterp.expectedPhigh(I));
        end

        figure; set(gcf,'Position',[86 925 1070 420]);
        subplot(1,2,1);
        % NOTE: data vals here come from dotsParse
        h(1) = errorbar(cohs, pRightHigh, pRightSEhigh, 'bo', 'MarkerFaceColor', 'b'); hold on;
        h(2) = errorbar(cohs, pRightLow, pRightSElow, 'bo', 'MarkerFaceColor', 'w');
        plot(cohs_fit,pRightHigh_model,'c-');
        plot(cohs_fit,pRightLow_model,'c--');
        xlabel('motion strength (%coh)');
        ylabel('proportion rightward choices');
        legend(h,'high bet','low bet','Location','Northwest');

        subplot(1,2,2); errorbar(cohs,pHigh,pHighSE,'ro');
        hold on; plot(cohs_fit,pHigh_model,'m--'); ylim([0 1]);
        xlabel('motion strength (%coh)');
        ylabel('proportion high bets');    
        
    case 2 %errfcn_DDM_2D_wConf
        
        ncohsteps = 99; % should be odd so there's a zero
        g_str = linspace(min(data.scoh),max(data.scoh),ncohsteps);
        nreps = 200;

        fitInterp = struct;
        fitInterp.scoh = repmat(g_str', nreps, 1);
        fitInterp.coherence = abs(fitInterp.scoh);

        % just fill these with ones, since error calc doesn't matter, and
        % they will be replaced with model values
        fitInterp.choice = ones(size(fitInterp.scoh));
        fitInterp.PDW = ones(size(fitInterp.scoh));
        fitInterp.correct = ones(size(fitInterp.scoh));
        fitInterp.RT = ones(size(fitInterp.scoh));

        options.plot = 0; options.feedback=0; fixed(:)=1;
        [~,~,fitInterp] = options.errfcn(X,X,fixed,fitInterp,options);
        
        cohs_fit = unique(fitInterp.scoh);
        n = nan(length(cohs_fit),1);
        pRightHigh_model = n;
        pRightLow_model = n;
        pHigh_model = n;
        RTmean_model = n;
        sigmaRT_model = n;

        % this parse + plot step could be streamlined, by making the above
        % scripts into functions that take data (or fitInterp) as argument
        for c = 1:length(cohs_fit)
            J = fitInterp.scoh==cohs_fit(c);

            nCor(c) = sum(J & fitInterp.correct); % use only correct trials for RT

            pRightHigh_model(c) = sum(J & fitInterp.choice==1 & fitInterp.PDW==1) / sum(J & fitInterp.PDW==1);
            pRightLow_model(c) = sum(J & fitInterp.choice==1 & fitInterp.PDW==0) / sum(J & fitInterp.PDW==0);
            pHigh_model(c) = sum(J & fitInterp.PDW==1) / sum(J);

            RTmean_model(c) = mean(fitInterp.RT(J & fitInterp.correct));
            sigmaRT_model(c) = std(fitInterp.RT(J & fitInterp.correct))/sqrt(nCor(c));
        end
        
        figure; set(gcf,'Position',[100 925 1470 420]);
        subplot(1,3,1);
        % NOTE: data vals here come from dotsParse
        h(1) = errorbar(cohs, pRightHigh, pRightSEhigh, 'bo', 'MarkerFaceColor', 'b'); hold on;
        h(2) = errorbar(cohs, pRightLow, pRightSElow, 'bo', 'MarkerFaceColor', 'w');
        plot(cohs_fit,pRightHigh_model,'c-');
        plot(cohs_fit,pRightLow_model,'c--');
        xlabel('motion strength (%coh)');
        ylabel('proportion rightward choices');
        legend(h,'high bet','low bet','Location','Northwest');

        subplot(1,3,2); errorbar(cohs,pHigh,pHighSE,'ro');
        hold on; plot(cohs_fit,pHigh_model,'m--'); ylim([0 1]);
        xlabel('motion strength (%coh)');
        ylabel('proportion high bets'); 
        
        subplot(1,3,3); errorbar(cohs,RTmean,RTse,'ro');
        hold on; plot(cohs_fit,RTmean_model,'m--');
        xlabel('motion strength (%coh)');
        ylabel('reaction time (s)'); 
end




%% testing out Gabe's code

% convert data struct to Gabe's format:

%   d = 
% 
%   1 x nStimulusConditions struct array with fields:
%     coherence   (1x1, signifies coherence for this stim cond.)
%     RTs         (nTrials x 1, Measured reaction time for each trial in this stim cond.)
%     choice      (nTrials x 1, Measured choice for each trial. 1 = right, 0 = left)
%     nTrials     (1x1, number of trials for this stimulus condition)
% 

u= unique(data.coherence);
for i=1:length(u)
  cohIdx{i}=find(data.coherence==u(i)); %get indeces into cell array
  data2(i).coherence = u(i); % one coherence value for this set of trials
  data2(i).RTs = data.RT(cohIdx{i}); % all of the trials for this coh
  data2(i).choice = data.choice(cohIdx{i}); % all choices for this coh
  data2(i).nTrials = length(cohIdx{i});
end




