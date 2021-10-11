function [err,fit] = dots3DMP_fit_2Dacc_err_sepbounds_noSim(param, guess, fixed, data, options)

% SJ 10-11-2021 no Monte Carlo simulation for fitting, it's redundant! just use model
% predictions directly

global call_num

% set parameters for this run based on guess and fixed params flag
param = getParam(param, guess, fixed);

if isfield(data,'PDW')
    data.conf = data.PDW;
end

mods   = unique(data.modality)';
cohs   = unique(data.coherence)'; 
hdgs   = unique(data.heading)';
% deltas = unique(data.delta)';
deltas = 0;

duration = 2000; % stimulus duration (ms)

kves  = param(1);
kvis  = param([2 3]);
sigmaVes = param(4);
sigmaVis = param([5 6]);
BVes     = abs(param(7)); % don't accept negative bound heights
BVis     = abs(param(8)); % fixed across cohs
BComb    = abs(param(9));
muTndVes = param(10);
muTndVis = param(11); % fixed across cohs
muTndComb = param(12);

if numel(param)==13, theta = param(13); end % only relevant for PDW

paramNames = {'kves','kvisLow','kvisHigh','sigmaVes','sigmaVisLow','sigmaVisHigh','BVes','BVis','BComb','muTndVes','muTndVis','muTndComb','theta'};

if options.conftask==2
    timeToConf = 350; % additional processing time for confidence
else
    timeToConf = 0;
end
duration = duration + timeToConf;

sdTnd = 60; % fixed SD
% Tnds = muTnd + randn(ntrials,1).*sdTnd; % fixed for all sims of given trial

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);

RVes.t = 0.001:0.001:duration/1000;
RVes.Bup = BVes;
RVes.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVes.lose_flag = 1;
RVes.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PVes =  images_dtb_2d(RVes);

RVis.t = 0.001:0.001:duration/1000;
RVis.Bup = BVis;
RVis.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVis.lose_flag = 1;
RVis.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PVis =  images_dtb_2d(RVis);

RComb.t = 0.001:0.001:duration/1000;
RComb.Bup = BComb;
RComb.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RComb.lose_flag = 1;
RComb.plotflag = 0; % 1 = plot, 2 = plot and export_fig

PComb =  images_dtb_2d(RComb);

%%

RTdata = data.RT*1000;
confdata = data.conf;

if options.conftask==1
    confdata = confdata*1000;
end

choiceD  = logical(data.choice-1);

pRight_model = nan(length(data.heading),1);
pHigh_model  = nan(length(data.heading),1);

n = nan(length(mods),length(cohs),length(deltas),length(hdgs));

meanRT_fit = n;     meanRT_data = n;    sigmaRT = n;
meanConf_fit = n;   meanConf_data = n;  sigmaConf = n;
nCor = n;

for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)

    for h = 1:length(hdgs)
        uh = abs(hdgs(h))==(hdgs(hdgs>=0));
        
        Jdata  = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        n(m,c,d,h) = sum(Jdata);

        % pRight_model stores the predictions made by the model for the
        % probability of a rightward choice on each trial condition in the
        % data
        % and then assigned to the actual data trials of each condition in pRight_model for
        % log likelihood calculation
        
        if m==1,     Ptemp = PVes; muTnd = muTndVes;
        elseif m==2, Ptemp = PVis; muTnd = muTndVis;
        elseif m==3, Ptemp = PComb; muTnd = muTndComb;
        end
        
        % SJ 10-11-2021 replaced redundant Monte Carlo simulation with computations
        % from method of images 
        % for choices, P.up/lo.p contains probabilities of bound hit
        % fo
        
        if hdgs(h)<0
            % if hdg is leftward, then pRight is p(incorrect response)
            pRight_model(Jdata) = Ptemp.lo.p(uh);  
        else
            pRight_model(Jdata) = Ptemp.up.p(uh);  
        end
        
        
        % RT/Conf Mean and SE from data, Conf Mean from model
        % use only correct trials (and 0
        N = 10000;
        conf_distr = Ptemp.up.pdf_t(uh,:) .* squeeze(Ptemp.up.distr_loser(uh,:,:))';
        conf_distr(conf_distr<0) = 0;
        try
            confSamples = randsample(Ptemp.logOddsCorrMap(:),N,true,conf_distr(:));
        catch
            keyboard
        end
        % use only correct trials + all zero heading trials (for RT and
        % SEP)
        
        usetrs_data  = data.correct | data.heading==0;
        
        nCor(m,c,d,h) = sum(Jdata & usetrs_data);
        
        if options.RTtask            
            meanRT_fit(m,c,d,h) = Ptemp.up.mean_t(uh)*1000+muTnd;
            
            meanRT_data(m,c,d,h) = mean(RTdata(Jdata & usetrs_data));
            sigmaRT(m,c,d,h) = std(RTdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
        end
        
        if options.conftask == 1 % sacc endpoint
            meanConf_fit(m,c,d,h) = 2*mean(logistic(confSamples))-1;
            
            meanConf_data(m,c,d,h) = mean(confdata(Jdata & usetrs_data));
            sigmaConf(m,c,d,h) = std(confdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
        elseif options.conftask == 2 % PDW
            pHigh_model(Jdata) = sum(confSamples > theta) / N;
        end
    end
    
end
end
end



%% calculate log likelihoods

% % likelihood of rightward choice on each trial, under binomial assumptions
% based on Palmer et al. 2005 coh version, but this ignores mod/coh in 3DMPtask
% Pr_model = 1 ./ (1 + exp(-2*k*B*sind(hdg)));

% kluge to avoid log(0) issues
pRight_model(pRight_model==0) = eps; 
pRight_model(pRight_model==1) = 1-eps;

% log likelihood of choice is summed log likelihood across all trials (log
% probability of observing choice given model)
LL_choice = (nansum(log(pRight_model(choiceD))) + nansum(log(1-pRight_model(~choiceD))));

% log likelihood of mean RTs, under Gaussian approximation
LL_RT = 0;
if options.RTtask
    L_RT = 1./(sigmaRT.*sqrt(2*pi)) .* exp(-(meanRT_fit - meanRT_data).^2 ./ (2*sigmaRT.^2));
    L_RT(L_RT==0) = eps; %min(L_RT(L_RT~=0));
    LL_RT = nansum(log(L_RT(:))); % sum over all conditions
end

% likelihood of mean confidence (for SEP), or log probabilities for PDW
LL_conf = 0; % initialize, if no conf task will stay as 0
if options.conftask == 1 % SEP
    L_conf = 1./(sigmaConf*sqrt(2*pi)) .* exp(-(meanConf_fit - meanConf_data).^2 ./ (2*sigmaConf.^2));
    L_conf(L_conf==0) = eps; %min(L_conf(L_conf~=0));
    LL_conf = nansum(log(L_conf(:))); % sum over all conditions
    
elseif options.conftask == 2 % PDW
    pHigh_model(pHigh_model==0) = eps; 
    pHigh_model(pHigh_model==1) = 1-eps;
    PDW = logical(data.conf);
    LL_conf = (nansum(log(pHigh_model(PDW))) + nansum(log(1-pHigh_model(~PDW))));   
end


err = -(LL_choice + LL_conf + LL_RT); % negate for minimization
% err = -(LL_choice + LL_conf);
% err = -LL_choice;

fit.pRight_model = pRight_model;
if options.conftask==1
    fit.meanConf_fit = meanConf_fit;
elseif options.conftask==2
    fit.pHigh_model = pHigh_model;
end
fit.meanRT_fit = meanRT_fit;



%% print progress report!
fprintf('\n\n\n****************************************\n');
fprintf('run %d\n', call_num);
fprintf('param values\n')
for i=1:length(param)
    fprintf('%s:\t%g\n',paramNames{i},param(i))
end
fprintf('\n');

fprintf('err: %f\n', err);
fprintf('split errs: Choice: %.2f, Conf: %.2f, RT: %.2f\n',LL_choice,LL_conf,LL_RT)
if options.ploterr && strcmp(options.fitMethod,'fms')
    figure(options.fh); hold on
    plot(call_num, err, '.','MarkerSize',12,'color','k');
    drawnow;
end
call_num = call_num + 1;


end


% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam ( param1 , guess , fixed )
  
  if all(fixed), param2 = param1; return; end

  param2(fixed==0) = param1;            %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess (initial point)
end
