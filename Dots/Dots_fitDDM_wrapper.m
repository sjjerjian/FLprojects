% Dots_fitDDM_wrapper.m

% Generalized wrapper script for Dots DDM fitting, formerly a part of
% Dots_offlineAnalysis.m

% requires a struct data with at minimum a variable for choice and one for
% signed coherence

% CF updated 12/2021

% % if running this standalone, load a data file, parse and plot it
clear all; close all;

% load Hanzo_data_fall2020.mat
load tempsim.mat

options.RTtask = 1;
options.conftask = 2;
cohs = unique(data.scoh);

% assign nearly-zero coh to exactly zero
data.coherence(data.coherence<1e-10) = 0;

% parse trial data into aggregated and other support vars
RTCorrOnly = 0;
if ~exist('parsedData','var')  % e.g., if simulation was run
    parsedData = Dots_parseData(data,options.conftask,options.RTtask,RTCorrOnly);
end

% optional [data will be plotted below regardless, along with the fits]
Dots_plot(parsedData,cohs,options.conftask,options.RTtask)



%% now the fitting itself


%****** first select which model to fit ********
% modelID=1; options.errfcn = @errfcn_DDM_1D_wConf;      % 1D DDM with threshold on log odds, usually for var dur [Kiani 09 (FP4)]
modelID=2; options.errfcn = @errfcn_DDM_2D_wConf_noMC; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (WolpertMOI)]

% options.errfcn = @errfcn_DDM_2D_wConf_noMC_signed; modelID=2; % as above, but with signed cohs
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
        if exist('origParams','var') % simulation
            k = origParams.k;
            B = origParams.B;
            theta = origParams.theta;
            alpha = origParams.alpha;
        else
            k = 0.6; % sensitivity parameter
            B = 15; % bound height
            theta = 0.8; % criterion (in log odds correct) for betting high
            alpha = 0.1; % base rate of low-bet choices
        end
        guess = [k B theta alpha];
        fixed = [0 0 0     0]; % can fix some params and fit the others, or fix all to hand-tune
        data.dur = round(data.duration*1000); % dur must be integer valued (in ms)
                             %^ same as RT for RT task
        
    case 2 %errfcn_DDM_2D_wConf
        
        if exist('origParams','var') % simulation
            k = origParams.k;
            B = origParams.B;
            sigma = origParams.sigma;
            theta = origParams.theta;
            alpha = origParams.alpha;
            Tnd = origParams.TndMean/1000; % convert to s
        else
            k = 20;
            B = 1;
            sigma = 0.05;
            theta = 2.0;
            alpha = 0; % base rate of low-bet choices
            Tnd = 0.3; % non-decision time (s)
        end
        
        guess = [k B theta alpha Tnd];
        fixed = [0 0 0     0     0];
        
%    case [others yet to be written: extrema, snapshot, SDT models?,...]

end

% ************************************
% set all fixed to 1 for hand-tuning:
fixed(:)=1;
% ************************************

% fit it!
[X, err_final, fit, parsedFit] = Dots_fitDDM(guess,fixed,data,options);


%% plot it!
Dots_plot_fit(parsedData,parsedFit,cohs,options.conftask,options.RTtask);


%% OLD:

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
        
% %         ncohsteps = 99; % should be odd so there's a zero
% %         g_str = linspace(min(data.scoh),max(data.scoh),ncohsteps); %
        % (this strategy works for 1D but in the current version of 2D
        % the model fits/predictions are generated by monte carlo
        % simulation and are too noisy to be worth repeating at
        % finely spaced headings. So for now hdgs will be the actual hdgs
        % from data although ncohsteps can increase in future iterations)
        g_str = unique(data.scoh)';
        
        nreps = 300;

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
        nCor = n;
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
            pHigh_model(c) = sum(J & fitInterp.PDW==1) / sum(J) - X(4); % X(4) is alpha!

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
