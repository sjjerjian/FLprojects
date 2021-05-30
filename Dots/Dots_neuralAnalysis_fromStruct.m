% NEW: load/plot data from Miguel struct


%%


%%


%%


%%


% clear
% close all

% load PrelimDataStruct
% load DataStruct_10-26-2020_12-18-2020

plotflag = 0;

extRaster = [50 150]; % extend raster this much before dots onset (1) and after dots offset (2)
latency = [80 120]; % cutoffs for spike rate [vs dots onset (1) and offset (2)]; latency(2) must be < extRaster(2)!
convKernel = fspecial('average', [1 40]); % N ms wide boxcar (acausal, centered)

% some manual plotting lims
psth_xlim1 = 400;
psth_xlim2 = [-300 100];

% for tuning functions
tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
tuning_vonMises_err = @(b,dir,FR) nansum((tuning_vonMises(b,dir)-FR).^2); % sum squared error

maxCoh = 0; % maximum coherence for inclusion in CP/confP and choice/conf-cond PSTH
% should be 0! unless testing code / sanity checks, or calcuating residuals


% quality control, etc

% check for missing spikes in mapping, use remap to recover if available,
% otherwise delete the channel [ONLY REQUIRED FOR ANALYSES THAT REQUIRE TUNING!]
deletedChan = nan(length(dataCell),6);
movedRemaptoMap = nan(length(dataCell),6);
for n = 1:length(dataCell)
    for c = length(dataCell{n}.Mapping.Channel):-1:1 % go backwards
        if isnan(dataCell{n}.Mapping.SpikeTimes{c})
            if isfield(dataCell{n},'Remap')                
%                 if isnan(dataCell{n}.Remap.SpikeTimes{c})
                    % no spikes in c for map or remap, delete this chan
                    dataCell{n}.Mapping.Channel(c)=[];                    
                    dataCell{n}.Mapping.SpikeTimes(c)=[];                    
                    dataCell{n}.Mapping.Selectivity(c)=[];
                    dataCell{n}.Remap.Channel(c)=[];                    
                    dataCell{n}.Remap.SpikeTimes(c)=[];                    
                    dataCell{n}.Remap.Selectivity(c)=[];
                    dataCell{n}.Exp.Channel(c)=[];                    
                    dataCell{n}.Exp.SpikeTimes(c)=[];                    
                    dataCell{n}.Exp.Selectivity(c)=[];
                    deletedChan(n,c) = 1;
%                 else % remap has spikes, replace mapping(c) with remap(c)               
                    % uncomment when remap spike times are corrected:
%                     dataCell{n}.Mapping.Channel(c) = dataCell{n}.Remap.Channel(c);                    
%                     dataCell{n}.Mapping.SpikeTimes(c) = dataCell{n}.Remap.SpikeTimes(c);                    
%                     dataCell{n}.Mapping.Selectivity(c) = dataCell{n}.Remap.Selectivity(c);
%                     movedRemaptoMap(n,c) = 1;
%                 end
            else % no remap, and no spikes in Mapping, so delete this chan
                dataCell{n}.Mapping.Channel(c)=[];                    
                dataCell{n}.Mapping.SpikeTimes(c)=[];                    
                dataCell{n}.Mapping.Selectivity(c)=[];
                dataCell{n}.Exp.Channel(c)=[];                    
                dataCell{n}.Exp.SpikeTimes(c)=[];                    
                dataCell{n}.Exp.Selectivity(c)=[];
                deletedChan(n,c) = 1;
            end
        end
    end
end


% check max trial len, nchans, ntrials (for preallocation)
for n = 1:length(dataCell)
    mStart = dataCell{n}.Mapping.openEvents.motionStart';
    mEnd = dataCell{n}.Mapping.openEvents.motionEnd';
    dur = round((mEnd-mStart)*1000);
    maxlen_map(n) = round(max(dur))+sum(extRaster);
    
    mStart = dataCell{n}.Exp.openEvents.motionStart';
    mEnd = dataCell{n}.Exp.openEvents.motionEnd';
    dur = round((mEnd-mStart)*1000);
    maxlen_exp(n) = round(max(dur))+sum(extRaster);
    
    nCh(n) = length(dataCell{n}.Mapping.Channel);
    nTrials_map(n) = length(dataCell{n}.Mapping.openEvents.motionStart);
    nTrials_exp(n) = length(dataCell{n}.Exp.openEvents.motionStart);
    
%     % this is to find any mismatches in channel vectors (other than NaNs)
%     if any(dataCell{n}.Mapping.Channel ~= dataCell{n}.Exp.Channel)
%         if ~isnan(dataCell{n}.Mapping.Channel(dataCell{n}.Mapping.Channel~=dataCell{n}.Exp.Channel)) && ...
%            ~isnan(dataCell{n}.Exp.Channel(dataCell{n}.Mapping.Channel~=dataCell{n}.Exp.Channel))     
%             keyboard
%         end
%     end
end
% this tells me that 1600 + sum(extRaster) is enough for exp,
% and 400 + sum(extRaster) is enough for mapping. set these manually and
% exclude trials that exceed them.
% the max nch and ntr will be set below.
rasterLen_exp = 1600;
rasterLen_map = 400;

% preallocate:
% raster: [sess,chan,trials,time]
raster_map = nan(size(dataCell,2),max(nCh),max(nTrials_map),rasterLen_map+sum(extRaster));
raster_dotsOn = nan(size(dataCell,2),max(nCh),max(nTrials_exp),rasterLen_exp+sum(extRaster));
raster_RT = raster_dotsOn;

% % aggregated, normalized psth [sess,chan,conditions(x4),time]:
% % psthNorm_dotsOn = nan(size(dataCell,2),max(nCh),4,rasterLen_exp+sum(extRaster));
% % psthNorm_RT = nan(size(dataCell,2),max(nCh),4,rasterLen_exp+sum(extRaster));
% ^ this wasn't easy to deal with, so decided to nindex by [sess*chan,time]
%   and therefore need one matrix for each condition (x2 for dotsOn and RT)
psthNorm_dotsOn_one = nan(sum(nCh),rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_two = nan(sum(nCh),rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_three = nan(sum(nCh),rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_four = nan(sum(nCh),rasterLen_exp+sum(extRaster));
psthNorm_RT_one = nan(sum(nCh),rasterLen_exp+sum(extRaster));
psthNorm_RT_two = nan(sum(nCh),rasterLen_exp+sum(extRaster));
psthNorm_RT_three = nan(sum(nCh),rasterLen_exp+sum(extRaster));
psthNorm_RT_four = nan(sum(nCh),rasterLen_exp+sum(extRaster));

% spike rates/counts: [sess,chan,trials]
spRate_map = nan(size(dataCell,2),max(nCh),max(nTrials_map));
spRate_exp = nan(size(dataCell,2),max(nCh),max(nTrials_exp));
spCount_exp = nan(size(dataCell,2),max(nCh),max(nTrials_exp));

% misc
prefDir_fit = nan(size(dataCell,2),max(nCh));
CP_high = nan(size(dataCell,2),max(nCh));
CP_low = nan(size(dataCell,2),max(nCh));
CP_all = nan(size(dataCell,2),max(nCh));
ConfP_pref = nan(size(dataCell,2),max(nCh));
ConfP_null = nan(size(dataCell,2),max(nCh));
ConfP_all = nan(size(dataCell,2),max(nCh));



% calc rasters, rates/counts, and psths for each session and channel
M=1;
for n = 1:length(dataCell)
    tic
    n
    %*****************************
    % TUNING (map)
    %*****************************
        
    mStart = dataCell{n}.Mapping.openEvents.motionStart';
    mEnd = dataCell{n}.Mapping.openEvents.motionEnd';
    dur = round((mEnd-mStart)*1000); % some variability (due to dropped frames???), but okay for now
    dir = dataCell{n}.Mapping.openEvents.direction;
    dirs = unique(dir);
    if any(movedRemaptoMap(n,:)==1)
        mStart_remap = dataCell{n}.Mapping.openEvents.motionStart';
        mEnd_remap = dataCell{n}.Mapping.openEvents.motionEnd';
        dur_remap = round((mEnd-mStart)*1000); % some variability (due to dropped frames???), but okay for now
        dir_remap = dataCell{n}.Mapping.openEvents.direction;
    end
    
    for c = 1:nCh(n)
        if movedRemaptoMap(n,c)==1 % if using spikes from remap, need to align sptimes with the corresponding events 
            SpikeTimes = dataCell{n}.Remap.SpikeTimes{c};
            startTrial = mStart_remap-extRaster(1)/1000;
            endTrial = mEnd_remap+extRaster(2)/1000;
            trialLen = endTrial - startTrial; % in seconds
            duration = dur_remap;
            nTr = length(dataCell{n}.Remap.openEvents.motionStart);
        else
            SpikeTimes = dataCell{n}.Mapping.SpikeTimes{c};
            startTrial = mStart-extRaster(1)/1000;
            endTrial = mEnd+extRaster(2)/1000;
            trialLen = endTrial - startTrial; % in seconds
            duration = dur;
            nTr = nTrials_map(n);
        end
        for t = 1:nTr
            if round(trialLen(t)*1000)<=size(raster_map,4)
                raster_map(n,c,t,1:round(trialLen(t)*1000)) = 0; % init bins to zero but keep nans after extRaster(2) cutoff
                spkInd = SpikeTimes >= startTrial(t) & SpikeTimes <= endTrial(t);
                spkTimes = SpikeTimes(spkInd) - startTrial(t);                
                if ~isempty(spkTimes)
                    x = 0:0.001:ceil(spkTimes(end)*1000)/1000;
                    raster_map(n,c,t,1:length(x)-1) = histcounts(spkTimes,x);
                    a = latency(1)+extRaster(1);
                    b = a + duration(t) - latency(1) + latency(2);
                    spRate_map(n,c,t) = sum(raster_map(n,c,t,a:b)) / ((b-a+1)/1000);
                end % else, will remain nans
            end % else, trial exceeds the max length specified in prealloc (will stay nans)
        end
    end

    % fit the tuning curves to get pref dir
    for c = 1:nCh(n)
        if ~isnan(dataCell{n}.Mapping.SpikeTimes{c})
            if movedRemaptoMap(n,c)==1 % if using spikes from remap, need to use dir from the corresponding events
                FRmean = calc_mean(squeeze(spRate_map(n,c,1:length(dir_remap))),dir_remap,dirs_remap);
                peak = dirs_remap(FRmean==max(FRmean));
                % vonMises params: b = [ampl, kappa, theta, offset]
                guess = [max(FRmean)-min(FRmean), 1.5, peak(1), min(FRmean)];
                [beta, fval] = fminsearch(@(x) tuning_vonMises_err(x,dir_remap,squeeze(spRate_map(n,c,1:length(dir_remap)))), guess);
            else
                FRmean = calc_mean(squeeze(spRate_map(n,c,1:length(dir))),dir,dirs);
                peak = dirs(FRmean==max(FRmean));
                % vonMises params: b = [ampl, kappa, theta, offset]
                guess = [max(FRmean)-min(FRmean), 1.5, peak(1), min(FRmean)];
% % %                 tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
% % %                 tuning_vonMises_err = @(b,dir,FR) nansum((tuning_vonMises(b,dir)-FR).^2); % sum squared error
                [beta, fval] = fminsearch(@(x) tuning_vonMises_err(x,dir,squeeze(spRate_map(n,c,1:length(dir)))), guess);
    %            % temp, to check fits by eye: (they're all good!)
    %             c
    % %             err = fval/length(dir);
    %             figure(1); clf; plot(dirs,FRmean,'bo-',dirs,tuning_vonMises(beta,dirs),'r-');
    %             pause
            end
            prefDir_fit(n,c) = beta(3);
        end
    end
  
    if plotflag
%     figure(n); set(gcf, 'Color', 'w', 'Position', [200 500 370*nCh(n) 770], 'PaperPositionMode', 'auto');
    figure(n); set(gcf, 'Color', 'w', 'Position', [1 185 370*nCh(n) 770], 'PaperPositionMode', 'auto');
    for c = 1:nCh(n)
        % simple tuning curve
        subplot(2,nCh(n),c);
        if movedRemaptoMap(n,c)==1 % if using spikes from remap, need to use dir from the corresponding events
            dirPlot = dir_remap; dirsPlot = dirs_remap;
        else
            dirPlot = dir; dirsPlot = dirs;
        end
        [FRmean, FRse, ~] = calc_mean(squeeze(spRate_map(n,c,1:length(dirPlot))),dirPlot,dirsPlot);
        errorbar(dirsPlot,FRmean,FRse,'ko-','MarkerFaceColor','k');
        set(gca,'xtick',dirsPlot,'xticklabel',{'0','','90','','180','','270',''},'tickdir','out')
        xlabel('direction (deg)'); ylabel('firing rate (sp/s)');
        title(sprintf('sess %d, chan %d',n,dataCell{n}.Mapping.Channel(c)));
        changeAxesFontSize(gca, 22, 22);
                
        % PSTH for pref/null, aligned to dots onset
        subplot(2,nCh(n),c+nCh(n));
        pref = dirsPlot(FRmean==max(FRmean)); null = dirsPlot(FRmean==min(FRmean));
        psth = calc_mean(squeeze(raster_map(n,c,:,:)), dirPlot, [pref;null])*1e3;
        psth = smoothRaster(psth, convKernel);
        tAxis = -extRaster(1) : size(squeeze(raster_map(n,c,:,:)),2) - extRaster(1) - 1;
        tAxis = tAxis + length(convKernel)/2+1; % offset for causal
        clr = {'g','r'};
        YLim = [min(min(psth))*0.9 max(max(psth))*1.1];
        XLim = [tAxis(1) 350];
        hold on;
        for p = 1:size(psth,1)
            plot(tAxis, psth(p,1:length(tAxis)),'Color', clr{p},'LineWidth',2);
        end
        set(gca, 'XLim', XLim, 'YLim', YLim, 'xtick', 0:100:400, 'xticklabel', {'0','','200','','400'}, 'TickDir','out'); box off;
        xlabel('Time from dots on (ms)');
        ylabel('Firing rate (spikes/s)');
        changeAxesFontSize(gca, 22, 22);
    end
    pause;
    end
    
    
    %*****************************
    % TASK (exp)
    %*****************************
        
    mStart = dataCell{n}.Exp.openEvents.motionStart';
        nTrials = length(mStart);
    mEnd = dataCell{n}.Exp.openEvents.motionEnd';
    RT = round((mEnd-mStart)*1000);
    dur=RT;
    dir = dataCell{n}.Exp.openEvents.direction;
        dirs = unique(dir);
        
    for c = 1:nCh(n)
        startTrial = mStart-extRaster(1)/1000;
        endTrial = mEnd+extRaster(2)/1000;
        trialLen = endTrial - startTrial; % in seconds
        for t = 1:nTrials_exp(n)
            % raster, align dotsOn (& spike counts/rates)
            if round(trialLen(t)*1000)<=size(raster_dotsOn,4)
                raster_dotsOn(n,c,t,1:round(trialLen(t)*1000)) = 0; % init bins to zero but keep nans after RT+extRaster(2)
                spkInd = dataCell{n}.Exp.SpikeTimes{c} >= startTrial(t) & dataCell{n}.Exp.SpikeTimes{c} <= endTrial(t);
                spkTimes = dataCell{n}.Exp.SpikeTimes{c}(spkInd) - startTrial(t);                
                if ~isempty(spkTimes)
                    x = 0:0.001:ceil(spkTimes(end)*1000)/1000;
                    raster_dotsOn(n,c,t,1:length(x)-1) = histcounts(spkTimes,x);
                    a = latency(1)+extRaster(1);
                    b = a + dur(t) - latency(1) + latency(2);
                    spRate_exp(n,c,t) = sum(raster_dotsOn(n,c,t,a:b)) / ((b-a+1)/1000);
                    spCount_exp(n,c,t) = sum(raster_dotsOn(n,c,t,a:b));
                end % else, spRate/Count will remain nans
            end % else, trial exceeds the max length specified in prealloc (will stay nans)
                        
            % raster, align RT (dotsOff)
            if round(trialLen(t)*1000)<=size(raster_RT,4)
                raster_RT(n,c,t,end-round(trialLen(t)*1000):end) = 0; % init bins to zero but keep nans before dotsOn-extRaster(1)
                spkTimes = flipud(endTrial(t)-dataCell{n}.Exp.SpikeTimes{c}(spkInd)); % flipud instead, makes hist easier
                if ~isempty(spkTimes)
                    x = 0:0.001:ceil(spkTimes(end)*1000)/1000;
                    raster_RT(n,c,t,end-length(x)+2:end) = fliplr(histcounts(spkTimes,x));
                end            
            end
            
            % identify pref choice targ
            I = dataCell{n}.Exp.openEvents.direction==0 & dataCell{n}.Exp.openEvents.coherence > 0.25; % high coh but some sess might exclude 0.512, 
            J = dataCell{n}.Exp.openEvents.direction==180 & dataCell{n}.Exp.openEvents.coherence > 0.25; % plus we get more s/n this way
            pref = double(mean(spRate_exp(n,c,I)) > mean(spRate_exp(n,c,J)));
        end
        
        % kluge to remove some anomalous high FRs (bimodal rate dist, obviously not real)
        if (n==1 && c==1) || (n==5 && c==1) % only seen it twice, so far
            K = spRate_exp(n,c,:)>110;
            raster_dotsOn(n,c,K,:) = nan;
            spCount_exp(n,c,K) = nan;
            spRate_exp(n,c,K) = nan;
        end
    end
    
    % CP, confP, and choice/conf-conditioned PSTHs
    J = abs(dataCell{n}.Exp.openEvents.coherence)<=maxCoh; % weak motion only, or else there's a confound! can also use residuals...

    % use rate or count for CP? you'd think count, because that's what the
    % rest of the brain has to work with, but then the choice-conditioned
    % distributions would not necessarily be equal under the null
    % hypothesis (on average they should be, with enough N, assuming equal
    % mean RTs for each choice -- but we can't assume that)
    for c = 1:nCh(n)        
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1;
            X = squeeze(spRate_exp(n,c,I)); % pref-high        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1;
            Y = squeeze(spRate_exp(n,c,I)); % null-high
        CP_high(n,c) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
                
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0;
            X = squeeze(spRate_exp(n,c,I)); % pref-low        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0;
            Y = squeeze(spRate_exp(n,c,I)); % null-low
        CP_low(n,c) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
        I = J & dataCell{n}.Exp.openEvents.choice==pref;
            X = squeeze(spRate_exp(n,c,I)); % pref-all        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref;
            Y = squeeze(spRate_exp(n,c,I)); % null-all
        CP_all(n,c) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1;
            X = squeeze(spRate_exp(n,c,I)); % pref-high        
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0;
            Y = squeeze(spRate_exp(n,c,I)); % pref-low
        ConfP_pref(n,c) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
% if n==17 && c==3
%     keyboard
%     figure; set(gcf, 'Color', 'w', 'Position', [800 100 450 300], 'PaperPositionMode', 'auto');
%     [Nx,~] = histcounts(X,0:5:70);
%     [Ny,E] = histcounts(Y,0:5:70);
%     bar(E(2:end)-2.5,[Nx;Ny]',.9); box('off');
%     set(gca,'YLim',[0 12],'ytick',0:2:12,'XLim',[0 70],'xtick',0:10:70,'TickDir','out');
%     legend('high bet','low bet'); legend('boxoff')
%     xlabel('Firing rate (sp/s)'); ylabel('Count');
%     changeAxesFontSize(gca, 16, 16);
%     export_fig('ConfP_example_dists', '-eps');
% end
        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1;
            X = squeeze(spRate_exp(n,c,I)); % null-high        
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0;
            Y = squeeze(spRate_exp(n,c,I)); % null-low
        ConfP_null(n,c) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
        
        I = J & dataCell{n}.Exp.openEvents.pdw==1;
            X = squeeze(spRate_exp(n,c,I)); % all-high        
        I = J & dataCell{n}.Exp.openEvents.pdw==0;
            Y = squeeze(spRate_exp(n,c,I)); % all-low
        ConfP_all(n,c) = rocN(X(~isnan(X)),Y(~isnan(Y)),100);
    end
    
    if plotflag
%         figure(n*10); set(gcf, 'Color', 'w', 'Position', [200 300 370*nCh(n) 430], 'PaperPositionMode', 'auto'); end
        figure(n*10); set(gcf, 'Color', 'w', 'Position', [1 525 370*nCh(n) 430], 'PaperPositionMode', 'auto');
    end
        % aligned dots-on
    for c = 1:nCh(n)
        clear psth;
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1; % pref-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(n,c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0; % pref-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(n,c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1; % null-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(n,c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0; % null-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_dotsOn(n,c,I,:)))*1e3, convKernel);

        rMax = max(max(psth(:,1:300)));
        psthNorm_dotsOn_one(M,:) = psth(1,:)/rMax;
        psthNorm_dotsOn_two(M,:) = psth(2,:)/rMax;
        psthNorm_dotsOn_three(M,:) = psth(3,:)/rMax;
        psthNorm_dotsOn_four(M,:) = psth(4,:)/rMax;

        if plotflag
            subplot(1,nCh(n),c);
            tAxis = -extRaster(1) : psth_xlim1;
            tAxis = tAxis + length(convKernel)/2+1; % offset for causal
            clr = {'g-', 'g--', 'r-', 'r--'};            
            YLim = [min(min(psth(:,1:length(tAxis))))*0.9 max(max(psth(:,1:length(tAxis))))*1.1];
            XLim = [tAxis(1) psth_xlim1];
            hold on;
            for p = 1:size(psth,1)
                plot(tAxis, psth(p,1:length(tAxis)), clr{p},'LineWidth',2);
            end
            set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
            xlabel('Time from dots on (ms)');
            ylabel('Firing rate (spikes/s)');
            changeAxesFontSize(gca, 22, 22);
            if c==nCh(n)
                legend('pref high','pref low','null high','null low','Location','Northeast'); legend('boxoff');
            end
        end
    end

        % aligned RT
    if plotflag
%         figure(n*100); set(gcf, 'Color', 'w', 'Position', [200 100 370*nCh(n) 430], 'PaperPositionMode', 'auto');
        figure(n*100); set(gcf, 'Color', 'w', 'Position', [1 20 370*nCh(n) 430], 'PaperPositionMode', 'auto');
    end
    for c = 1:nCh(n)
        clear psth
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==1; % pref-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_RT(n,c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==pref & dataCell{n}.Exp.openEvents.pdw==0; % pref-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_RT(n,c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==1; % null-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_RT(n,c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice~=pref & dataCell{n}.Exp.openEvents.pdw==0; % null-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_RT(n,c,I,:)))*1e3, convKernel);
        
        psthNorm_RT_one(M,:) = psth(1,:)/rMax;
        psthNorm_RT_two(M,:) = psth(2,:)/rMax;
        psthNorm_RT_three(M,:) = psth(3,:)/rMax;
        psthNorm_RT_four(M,:) = psth(4,:)/rMax;
        M=M+1

        if plotflag
            subplot(1,nCh(n),c);
            tAxis = -(size(raster_RT,4)-extRaster(2)-1) : extRaster(2);
            tAxis = tAxis + length(convKernel)/2+1; % offset for causal
            XLim = psth_xlim2;
            YLim = [min(min(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*0.9 max(max(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*1.1];
            hold on;
            for p = 1:size(psth,1)
                plot(tAxis, psth(p,:), clr{p},'LineWidth',2);
            end
            set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
            xlabel('Time from saccade (ms)');
            changeAxesFontSize(gca, 22, 22);
        end
    end
    
    toc
    if plotflag; pause; end
end

%% plot average choice/conf-conditioned PSTH

    % aligned dotsOn
clear psth
psth(1,:) = nanmean(psthNorm_dotsOn_one);
psth(2,:) = nanmean(psthNorm_dotsOn_two);
psth(3,:) = nanmean(psthNorm_dotsOn_three);
psth(4,:) = nanmean(psthNorm_dotsOn_four);

figure(5000); set(gcf, 'Color', 'w', 'Position', [1 220 700 450], 'PaperPositionMode', 'auto');
tAxis = -extRaster(1) : psth_xlim1;
tAxis = tAxis + length(convKernel)/2+1; % offset for causal
clr = {'b-', 'b--', 'r-', 'r--'};            
YLim = [min(min(psth(:,1:length(tAxis))))*0.9 max(max(psth(:,1:length(tAxis))))*1.1];
XLim = [tAxis(1) psth_xlim1];
hold on;
for p = 1:size(psth,1)
    plot(tAxis, psth(p,1:length(tAxis)), clr{p},'LineWidth',2);
end
set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
xlabel('Time from dots on (ms)');
ylabel('Firing rate (normalized)');
changeAxesFontSize(gca, 22, 22);
h = legend('pref choice, high bet','pref choice, low bet','null choice, high bet','null choice, low bet','Location','Northeast'); legend('boxoff');
% export_fig('psthNorm_dotsOn_coh0', '-eps');
% export_fig('psthNorm_dotsOn_coh32', '-eps');

    % aligned RT
clear psth
psth(1,:) = nanmean(psthNorm_RT_one);
psth(2,:) = nanmean(psthNorm_RT_two);
psth(3,:) = nanmean(psthNorm_RT_three);
psth(4,:) = nanmean(psthNorm_RT_four);

figure(5001); set(gcf, 'Color', 'w', 'Position', [350 220 375 450], 'PaperPositionMode', 'auto');
% subplot(1,2,2);
tAxis = -(size(raster_RT,4)-extRaster(2)-1) : extRaster(2);
tAxis = tAxis + length(convKernel)/2+1; % offset for causal
XLim = [-200 20];
YLim = [min(min(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*0.9 max(max(psth(:,end-(psth_xlim2(2)-psth_xlim2(1)):end)))*1.1];
hold on;
for p = 1:size(psth,1)
    plot(tAxis, psth(p,:), clr{p},'LineWidth',2);
end
set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
xlabel('Time from saccade (ms)');
changeAxesFontSize(gca, 22, 22);
% export_fig('psthNorm_RT_coh0', '-eps');
% export_fig('psthNorm_RT_coh0_extended', '-eps');
% export_fig('psthNorm_RT_coh32', '-eps');
% export_fig('psthNorm_RT_coh32_extended', '-eps');


%% CP and ConfP, vs each other and as a func of prefDir

nanmean(CP_high(:))
nanmean(CP_low(:))
nanmean(CP_all(:))
nanmean(ConfP_pref(:))
nanmean(ConfP_null(:))
nanmean(ConfP_all(:))


%% choose

% CP = CP_high(~isnan(CP_high));
% or
% CP = CP_low(~isnan(CP_low));
% or
CP = CP_all(~isnan(CP_all));

ConfP = ConfP_pref(~isnan(ConfP_pref));
% or
% ConfP = ConfP_null(~isnan(ConfP_null));
% or
% ConfP = ConfP_all(~isnan(ConfP_all));

CPse = std(CP)/sqrt(length(CP));
[~,P] = ttest(CP,0.5)

ConfPse = std(ConfP)/sqrt(length(ConfP))
[~,P] = ttest(ConfP,0.5)



%%
figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
plot([0 1],[0 1],'k--',[0 1],[0.5 0.5],'k:',[0.5 0.5],[0 1],'k:');
hold on; plot(CP,ConfP,'ko','MarkerSize',11,'MarkerFaceColor',[1 1 1]); axis square;
% set(gca,'XLim',[0.17 0.83],'YLim',[0.17 0.83],'xtick',0.2:0.1:0.8,'ytick',0.2:0.1:0.8,'TickDir','out')
set(gca,'XLim',[0.1 0.9],'YLim',[0.1 0.9],'xtick',0.1:0.2:0.9,'ytick',0.1:0.2:0.9,'TickDir','out')
xlabel('Choice probability'); ylabel('Confidence probability');
changeAxesFontSize(gca, 22, 22);
% export_fig('CPvsConfP_coh0', '-eps');

figure; set(gcf, 'Color', 'w', 'Position', [800 200 450 80], 'PaperPositionMode', 'auto');
[histCP,E] = histcounts(CP,0:0.05:1);
bar(E(2:end)-0.025,histCP); box('off');
set(gca,'XLim',[0.1 0.9],'xtick',0.1:0.2:0.9,'xticklabel',[],'TickDir','out');
changeAxesFontSize(gca, 14, 14);
% export_fig('CPbar_coh0', '-eps');

figure; set(gcf, 'Color', 'w', 'Position', [800 100 80 450], 'PaperPositionMode', 'auto');
[histCP,E] = histcounts(ConfP,0:0.05:1);
barh(E(2:end)-0.025,histCP); box('off');
set(gca,'YLim',[0.1 0.9],'ytick',0.1:0.2:0.9,'yticklabel',[],'TickDir','out');
changeAxesFontSize(gca, 14, 14);
% export_fig('ConfPbar_coh0', '-eps');

[r,p] = corrcoef(CP,ConfP)


%%

prefDir = prefDir_fit(~isnan(prefDir_fit));

figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
plot(prefDir,CP,'ko','MarkerSize',11,'MarkerFaceColor',[1 1 1]);
set(gca,'xlim',[0 360],'xtick',0:45:360,'tickDir','out')
xlabel('prefDir'); ylabel('CP');

figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
plot(prefDir,ConfP,'ko','MarkerSize',11,'MarkerFaceColor',[1 1 1]);
set(gca,'xlim',[0 360],'xtick',0:45:360,'tickDir','out')
xlabel('prefDir'); ylabel('ConfP');
% ^too noisy to see anything here?

I = cosd(prefDir)>sqrt(2)/2 | cosd(prefDir)<-sqrt(2)/2; % within +/- 45 deg of 0 or 180
mean(CP(I))
mean(CP(~I))
%^ still nothing!


%%
% running mean
windowWidth = 90;

figure; set(gcf, 'Color', 'w', 'Position', [400 200 450 450], 'PaperPositionMode', 'auto');
[x_,y_,~] = running_mean(prefDir, CP, windowWidth, 'window_width'); h(1) = plot(x_,y_,'k-','LineWidth',2);
[x_,y_,~] = running_mean(prefDir, ConfP, windowWidth, 'window_width'); hold on; h(2) = plot(x_,y_,'r-','LineWidth',2);
plot([0 360],[0.5 0.5],'k--');
set(gca,'xlim',[0 360],'xtick',0:45:360,'ylim',[0.45 0.6], 'ytick',0.45:0.05:0.6,'tickDir','out')
xlabel('Preferred direction (deg)'); ylabel('Choice/Confidence Probability');
legend(h(1:2),'CP','ConfP','location','northwest'); legend('boxoff');
changeAxesFontSize(gca, 16, 16);

% export_fig('CP-v-prefdir_runningMean', '-eps');

