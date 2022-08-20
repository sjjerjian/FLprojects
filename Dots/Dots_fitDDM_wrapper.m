% Dots_fitDDM_wrapper.m

% Generalized wrapper script for Dots DDM fitting, formerly a part of
% Dots_offlineAnalysis.m

% requires a struct data with at minimum a variable for choice and one for
% signed coherence

% CF updated 12/2021, again in 07/2022

% % if running this standalone, load a data file, parse and plot it
clear all; close all;

% load Hanzo_data_fall2020.mat
load tempsim.mat % temp, for (pre)param recovery
% load doubtconf.mat 

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

% % optional [data will be plotted below regardless, along with the fits]
% Dots_plot(parsedData,cohs,options.conftask,options.RTtask)


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
% options.fitMethod = 'fmc'; % fmincon
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% % options.fitMethod = 'bads'; % not implemented yet, see dots3DMP for work-in-progress

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
            theta = 2.6; % criterion (in log odds correct) for betting high
            alpha = 0.2; % base rate of low-bet choices
            Tnd = 1.5; % non-decision time (s)
        end
        guess = [k B theta alpha Tnd];
        fixed = [0 0 0     0     0  ]; % can fix some params and fit the others, or fix all to hand-tune
        
    case 2 %errfcn_DDM_2D_wConf
        
        if exist('origParams','var') % i.e., from simulation
            k = origParams.k;
            B = origParams.B;
            theta = origParams.theta;
            alpha = origParams.alpha;
            Tnd = origParams.TndMean/1000; % convert to s
        else
%             k = 20;
%             B = 1;
%             theta = 2.0;
%             alpha = 0;
%             Tnd = 0.3;
            
            % these are attempts at fitting DoubtConf
            k= 3;
        	B= 1.5;
            theta= 1.5;
            alpha= 0.1;
            Tnd= 0.8;
%             best err: 5686.443107

        end
        guess = [k B theta alpha Tnd];
        fixed = [0 0 0     0     0];
        
%    case [others yet to be written: extrema, snapshot, SDT models?,...]

end

% ************************************
% set all fixed to 1 for hand-tuning or "pre-param-recovery"
fixed(:)=1;

% Randomize guess, for true param recovery
% guess = guess.*(rand(1,length(guess))+0.5);
% ************************************

% fit it!
[X, err_final, fit, parsedFit] = Dots_fitDDM(guess,fixed,data,options);


% plot it!
Dots_plot_fit(parsedData,parsedFit,cohs,options.conftask,options.RTtask);




