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

if options.dummyRun
    deltas = unique(data.delta)';
else
    deltas = 0;
end

duration = 2000; % stimulus duration (ms)

kves  = param(1);
kvis  = param([2 3]);
BVes     = abs(param(4)); % don't accept negative bound heights
BVis     = abs(param(5)); % fixed across cohs
BComb    = abs(param(6));
muTnd    = param(7);
ttc      = param(8);

try theta = param(9); catch; end % only relevant for PDW

paramNames = {'kves','kvisLo','kvisHi','BVes','BVis','BComb','muTnd','T2Conf','theta'};

duration = duration + ttc;

sdTnd = 60; % fixed SD

% assume the mapping for confidence is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k is the
% mean
% still need separate logOddsMaps if bound heights are different
k = mean([kves kvis]);

RVes.t = 0.001:0.001:duration/1000;
RVes.Bup = BVes;
RVes.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVes.lose_flag = 1;
RVes.plotflag = 0; % 1 = plot, 2 = plot and export_fig
PVesConf =  images_dtb_2d(RVes);
VesLogOdds = PVesConf.logOddsCorrMap;

RVis.t = 0.001:0.001:duration/1000;
RVis.Bup = BVis;
RVis.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVis.lose_flag = 1;
RVis.plotflag = 0; % 1 = plot, 2 = plot and export_fig
PVisConf =  images_dtb_2d(RVis);
VisLogOdds = PVisConf.logOddsCorrMap;

RComb.t = 0.001:0.001:duration/1000;
RComb.Bup = BComb;
RComb.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RComb.lose_flag = 1;
RComb.plotflag = 0; % 1 = plot, 2 = plot and export_fig
PCombConf =  images_dtb_2d(RComb);
CombLogOdds = PCombConf.logOddsCorrMap;

% now compute the P we need for model choices and RTs (modality-specific)

RVes.t = 0.001:0.001:duration/1000;
RVes.Bup = BVes;
RVes.drift = kves * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
RVes.lose_flag = 1;
RVes.plotflag = 0; % 1 = plot, 2 = plot and export_fig
PVes =  images_dtb_2d(RVes);

for c = 1:length(cohs)
    RVis.t = 0.001:0.001:duration/1000;
    RVis.Bup = BVis;
    RVis.drift = kvis(c) * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
    RVis.lose_flag = 1;
    RVis.plotflag = 0; % 1 = plot, 2 = plot and export_fig
    PVis{c} =  images_dtb_2d(RVis);
    
    kcombsq = kves^2 + kvis(c)^2;
    
    wVes = sqrt( kves^2 / kcombsq );
    wVis = sqrt( kvis(c)^2 / kcombsq );
    
    for d = 1:length(deltas)
        % positive delta defined as ves to the left, vis to the right
        muVes = kves    * sind(hdgs(hdgs>=0)-deltas(d)/2);
        muVis = kvis(c) * sind(hdgs(hdgs>=0)+deltas(d)/2);
        
        RComb.t = 0.001:0.001:duration/1000;
        RComb.Bup = BComb;
        
%         kcomb(c) = sqrt(kcombsq);
%         RComb.drift = kcomb(c) * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
        
        Rcomb.drift = wVes.*muVes + wVis.*muVis;

        RComb.lose_flag = 1;
        RComb.plotflag = 0; % 1 = plot, 2 = plot and export_fig
        PComb{c,d} =  images_dtb_2d(RComb);
    end
end



%%

RTdata = data.RT*1000;
confdata = data.conf;

if options.conftask==1
    confdata = confdata*1000;
end

usetrs_data  = data.correct | data.heading==0; % use only correct (and 0 hdg) trials for RT/SEP fits
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
        
        if sum(Jdata)==0, continue, end
        
        n(m,c,d,h) = sum(Jdata);

        % pRight_model stores the predictions made by the model for the
        % probability of a rightward choice on each trial condition in the
        % data
        % and then assigned to the actual data trials of each condition in pRight_model for
        % log likelihood calculation
        
        if m==1,     Ptemp = PVes;  logOddsMap = VesLogOdds; %k = kves; B = BVes;
        elseif m==2, Ptemp = PVis{c};  logOddsMap = VisLogOdds; %k = kvis(c); B = BVis;
        elseif m==3, Ptemp = PComb{c,d}; logOddsMap = CombLogOdds; %k = kcomb(c); B = BComb;
        end
        
        % SJ 10-11-2021 replaced redundant Monte Carlo simulation with computations
        % from method of images 
        % for choices, P.up/lo.p contains probabilities of bound hit
        
        if hdgs(h)<0
            % if hdg is leftward, then pRight is p(incorrect response)
            pRight_model(Jdata) = Ptemp.lo.p(uh);  
        else
            pRight_model(Jdata) = Ptemp.up.p(uh);  
        end
        pRight_fit(m,c,d,h) = mean(pRight_model(Jdata));

        
        % Palmer 2005 method
%         pRight_model(Jdata) = 1 ./ (1 + exp(-2*k*B*sind(hdgs(h))));
        
%         if hdgs(h)==0,keyboard,end
        % RT/Conf Mean and SE from data, Conf Mean from model
%         N = 10000;
%         conf_distr = Ptemp.up.pdf_t(uh,:) .* squeeze(Ptemp.up.distr_loser(uh,:,:))';
        conf_distr = Ptemp.up.pdf_t(uh,:) .* circshift(squeeze(Ptemp.up.distr_loser(uh,:,:))',-ttc,2); % shift distr_loser forward in time by ttc
        conf_distr = conf_distr ./ sum(conf_distr(:)); % normalize to posterior
        conf_distr(conf_distr<0) = 0;
        %             confSamples = randsample(logOddsMap(:),N,true,conf_distr(:));

        nCor(m,c,d,h) = sum(Jdata & usetrs_data);
        
        % CONF
        if options.conftask == 1 % sacc endpoint
%             meanConf_fit(m,c,d,h) = mean(2*logistic(confSamples)-1);
            meanConf_fit(m,c,d,h) = nansum(conf_distr(:) .* (2*logistic(logOddsMap(:))-1));
            
            meanConf_data(m,c,d,h) = mean(confdata(Jdata & usetrs_data));
            sigmaConf(m,c,d,h) = std(confdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
            
        elseif options.conftask == 2 % PDW
%             pHigh_model(Jdata) = sum(confSamples > theta) / N;

            % integrate distribution where logOdds>theta
            temp  = conf_distr(:);
            isLow = logOddsMap(:)<=theta;
            pHigh_model(Jdata) = 1-nansum(temp(isLow));
            
            % store by condition for ease of plotting later
            pHigh_fit(m,c,d,h) = mean(pHigh_model(Jdata)); 
        end
        
        % RT
        if options.RTtask            
            meanRT_fit(m,c,d,h) = Ptemp.up.mean_t(uh)*1000+muTnd;
            
            meanRT_data(m,c,d,h) = mean(RTdata(Jdata & usetrs_data));
            sigmaRT(m,c,d,h) = std(RTdata(Jdata & usetrs_data)) / sqrt(nCor(m,c,d,h));
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


% assign same conditions for fit data, will be useful for plotting later!
fit.heading = data.heading;
fit.modality = data.modality;
fit.coherence = data.coherence;
fit.delta = data.delta;

fit.pRight = pRight_fit;
if options.conftask==1
    fit.confMean = meanConf_fit;
elseif options.conftask==2
    fit.confMean = pHigh_fit;
end
fit.RTmean = meanRT_fit/1000;


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
    if call_num==1,clf; end
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
