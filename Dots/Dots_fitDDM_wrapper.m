% Dots_fitDDM_wrapper.m

% Generalized wrapper script for Dots DDM fitting, formerly a part of
% Dots_offlineAnalysis.m

% requires a struct data with at minimum a variable for choice and one for
% signed coherence

% CF updated 12/2021, again in 07/2022

clear all; close all;

% datafile = 'tempsim.mat'; % simulated data, e.g. for (pre)param recovery
% datafile = 'doubtconf.mat'; % 'Doubt' dataset, Marton et al.
datafile = 'Hanzo_data_fall2020_withFall2021.mat';
% datafile = 'Genji_data';


load(datafile)

if ~isfield(data,'conf')
    data.conf = data.PDW;
end




%% some bookkeeping, then parse data, and plot if desired

if ~exist('allowNonHB','var'); allowNonHB=0; end

if allowNonHB==0
% ignore unabsorbed probability, ie no need for weighted sum with max_dur
% in RT calculation, e.g. when sim excluded trials that failed to hit bound
    options.ignoreUnabs = 1;
else
    options.ignoreUnabs = 0;    
end

options.RTtask = 1;
options.conftask = 2; % 1=continuous/rating, 2=PDW
cohs = unique(data.scoh);

% parse trial data into aggregated and other support vars
RTCorrOnly = 0;
if ~exist('parsedData','var')  % e.g., if simulation was run
    parsedData = Dots_parseData(data,options.conftask,options.RTtask,RTCorrOnly);
end

if ~exist('data.PDW_preAlpha','var')
    data.PDW_preAlpha = data.PDW;
end

% % ********
% % optional [data will be plotted below regardless, along with the fits]
forTalk = 0;
Dots_plot(parsedData,cohs,options.conftask,options.RTtask,0,forTalk)
% % ********



%% now the fitting itself


%****** first select which model to fit ********
% modelID=1; options.errfcn = @errfcn_DDM_1D_wConf;    % 1D DDM with threshold on log odds, usually for var dur [Kiani 09 (uses Kiani's FP4 (Chang-Cooper method))]

modelID=2; options.errfcn = @errfcn_DDM_2D_wConf_noMC; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (uses Wolpert's images_dtb_2d (method of images, from Moreno-Bote 2010))]
% modelID=2; options.errfcn = @errfcn_DDM_2D_wConf_noMC_signed; % as above, but with signed cohs
%***********************************************

if exist('confModel','var') % i.e., from simulation
    options.confModel = confModel;
else
    options.confModel  = 'evidence+time'; % default, Kiani14: theta is criterion in log odds correct for betting high
%     options.confModel  = 'evidence_only';
                    % here, theta is treated as a flat criterion on evidence,
                    % and is therefore in DV units (distance from bound, so
                    % that the positive association with conf is preserved)
%     options.confModel  = 'time_only';
                    % here, theta is simply a deadline: if past it, bet low
end

options.feedback = 1; % 1 = text output to cmd window, 2 = that, plus plot LL across runs
options.plot = 0; % plot the marginal PDFs, logOddsCorr map, and high/low bet regions (only makes sense for fixed(:)=1)

% choose optimization method
options.fitMethod = 'fms'; % fminsearch
% options.fitMethod = 'fmc'; % fmincon

% the next three are from global optim toolbox, but are very slow and not clearly better
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';

% Bayesian adaptive directed search, from Luigi Acerbi
% options.fitMethod = 'bads'; % not implemented yet, see dots3DMP for work-in-progress


% params: 
switch modelID
    case 1 %errfcn_DDM_1D_wConf
        
        % initial guess (or hand-tuned params)
        if exist('origParams','var') % i.e., from simulation
            k = origParams.k;
            B = origParams.B;
            theta = origParams.theta;
            alpha = origParams.alpha;
            Tnd = origParams.TndMean/1000; % convert to s
        else
            % these are attempts at fitting DoubtConf
            k = 0.3; % sensitivity parameter
            B = 40; % bound height
            theta = 2.6; % criterion for betting high
            alpha = 0.2; % base rate of low-bet choices
            Tnd = 1.5; % non-decision time (s)
        end
        guess = [k B theta alpha Tnd];
        fixed = [0 0 0     0     0  ]; % can fix some params and fit the others, or fix all to hand-tune
        
    case 2 %errfcn_DDM_2D_wConf
        
        % initial guess (or hand-tuned params)
        if exist('origParams','var') % i.e., from simulation
            k = origParams.k;
            B = origParams.B;
            theta = origParams.theta;
            alpha = origParams.alpha;
            try % sometimes want separate non-decision time for R and L
                TndR = origParams.TndRMean/1000; % convert to s
                TndL = origParams.TndLMean/1000;
            catch
                TndR = origParams.TndMean/1000;
                TndL = origParams.TndMean/1000;
            end
        else

            
            switch datafile(1:5)

                case 'doubt' % 'Doubt' dataset, Marton et al.
                    k= 4;
                    B= 2.15;
                    theta= 2.5;
%                     theta= 3.5;
                    alpha= 0.12;
                    TndR= 0.58;
                    TndL= 0.58;
                case 'Hanzo'
                    k= 26;
                    B= .52;
                    theta= 0.99;
                    alpha= 0.14;
                    TndR= 0.30;
                    TndL= 0.37;
                case 'Genji'
%                     % GENJI: 'fit' choice
%                     k= 10;
%                     B= .76;
%                     theta= 1.22;
%                     alpha= 0;
%                     TndR= 0.22;
%                     TndL= 0.22;                

%                      % GENJI: 'fit' RT
%                     k= 14;
%                     B= .74;
%                     theta= 1.35;
%                     alpha= 0;
%                     TndR= 0.27;
%                     TndL= 0.26;

                    % GENJI: 'fit' Conf
                    k= 19;
                    B= .74;
                    theta= 1.35;
                    alpha= 0;
                    TndR= 0.28;
                    TndL= 0.27;
    
            end

        end
        
        guess = [k B theta alpha TndR TndL];
        fixed = [0 0 0     0     0    0];
        
%    case [others yet to be written: extrema, snapshot, SDT models?,...]

end

% ************************************
% set all fixed to 1 for hand-tuning or "pre-param-recovery"
fixed(:)=1;

% % temp: directed perturbations to see if effect on LL is sensible
% guess(3) = 0.7; % affects LL for choice_low more than _high, although not as big a diff as might expect (also affects conf and RT, both curves)
% guess(4) = 0; % affects LL for conf only
% guess(5) = 0.6; guess(6) = 0.6; % affects LL for RT only


% Randomize guess, for true param recovery
% guess = guess.*(rand(1,length(guess))+0.5);
% ************************************

% fit it!
[X, err_final, fit, parsedFit] = Dots_fitDDM(guess,fixed,data,options);


% plot it!
forTalk = 0;
Dots_plot_fit(parsedData,parsedFit,cohs,options.conftask,options.RTtask,forTalk);


