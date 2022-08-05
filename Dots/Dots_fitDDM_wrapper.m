% Dots_fitDDM_wrapper.m

% Generalized wrapper script for Dots DDM fitting, formerly a part of
% Dots_offlineAnalysis.m

% requires a struct data with at minimum a variable for choice and one for
% signed coherence

% CF updated 12/2021, again in 07/2022

% % if running this standalone, load a data file, parse and plot it
clear all; close all;

% load Hanzo_data_fall2020.mat
load tempsim.mat % temp, for param recovery

options.RTtask = 1;
options.conftask = 2; % 1=continuous, 2=PDW
cohs = unique(data.scoh);

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
% modelID=2; options.errfcn = @errfcn_DDM_2D_wConf_noMC; % 2D DDM aka anticorrelated race, for RT+conf [Kiani 14 / van den Berg 16 (WolpertMOI)]

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
        if exist('origParams','var') % i.e., from simulation
            k = origParams.k;
            B = origParams.B;
            theta = origParams.theta;
            alpha = origParams.alpha;
            Tnd = origParams.TndMean/1000; % convert to s
        else
            k = 0.4; % sensitivity parameter
            B = 30; % bound height
            theta = 0.8; % criterion (in log odds correct) for betting high
            alpha = 0.1; % base rate of low-bet choices
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
            k = 20;
            B = 1;
            theta = 2.0;
            alpha = 0;
            Tnd = 0.3;
        end
        guess = [k B theta alpha Tnd];
        fixed = [0 0 0     0     0];
        
%    case [others yet to be written: extrema, snapshot, SDT models?,...]

end

% ************************************
% set all fixed to 1 for hand-tuning
% or "pre-parameter-recovery"
fixed(:)=1;
% ************************************

% fit it!
[X, err_final, fit, parsedFit] = Dots_fitDDM(guess,fixed,data,options);


% plot it!
Dots_plot_fit(parsedData,parsedFit,cohs,options.conftask,options.RTtask);




