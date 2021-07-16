% NEW: load/plot data from Miguel struct

% 5-2-21 : simpler version, spike counts/rates only

clear
close all

load DataStruct_10-26-2020_12-18-2020

plotflag = 0;

extRaster = [50 150]; % extend raster this much before dots onset (1) and after dots offset (2)
latency = [80 120]; % cutoffs for spike rate [vs dots onset (1) and offset (2)]; latency(2) must be < extRaster(2)!
convKernel = fspecial('average', [1 40]); % kernel for smoothing PSTHs (N ms wide boxcar, acausal, centered)

% for fitting tuning functions
tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
tuning_vonMises_err = @(b,dir,FR) nansum((tuning_vonMises(b,dir)-FR).^2); % sum squared error


%% quality control

% check for missing spikes in mapping, use remap to recover if available,
% otherwise delete the channel
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
% exclude trials that exceed them (monkey failed to fixate or some other
% glitch). the max nch and ntr will be set below.
rasterLen_exp = 1600;
rasterLen_map = 400;

%% preallocate:

% raster: [session,channel,trials,time]
raster_map = nan(size(dataCell,2),max(nCh),max(nTrials_map),rasterLen_map+sum(extRaster));
% % raster_dotsOn = nan(size(dataCell,2),max(nCh),max(nTrials_exp),rasterLen_exp+sum(extRaster));
% % raster_RT = raster_dotsOn;
  % ^ for expt data, may need later

% spike rates/counts: [session,channel,trials]
spRate_map = nan(size(dataCell,2),max(nCh),max(nTrials_map));
% % spRate_exp = nan(size(dataCell,2),max(nCh),max(nTrials_exp));
% % spCount_exp = nan(size(dataCell,2),max(nCh),max(nTrials_exp));
  % ^ for expt data, may need later

% misc
prefDir_fit = nan(size(dataCell,2),max(nCh));


%% calc rasters and spike rates/counts for each session and channel

M=1;
for n = 1:length(dataCell)
    tic
    n
    %*****************************
    % TUNING (map)
    %*****************************
        
    mStart = dataCell{n}.Mapping.openEvents.motionStart';
    mEnd = dataCell{n}.Mapping.openEvents.motionEnd';
    dur = round((mEnd-mStart)*1000); % some variability (due to dropped frames? more likely network lag), but okay for now
    dir = dataCell{n}.Mapping.openEvents.direction;
    dirs = unique(dir);
    
    for c = 1:nCh(n)
        SpikeTimes = dataCell{n}.Mapping.SpikeTimes{c};
        startTrial = mStart-extRaster(1)/1000;
        endTrial = mEnd+extRaster(2)/1000;
        trialLen = endTrial - startTrial; % in seconds
        duration = dur;
        nTr = nTrials_map(n);
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
            FRmean = calc_mean(squeeze(spRate_map(n,c,1:length(dir))),dir,dirs);
            peak = dirs(FRmean==max(FRmean));
                % vonMises params: b = [ampl, kappa, theta, offset]
            guess = [max(FRmean)-min(FRmean), 1.5, peak(1), min(FRmean)];
            [beta, fval] = fminsearch(@(x) tuning_vonMises_err(x,dir,squeeze(spRate_map(n,c,1:length(dir)))), guess);
%            % uncomment to check fits by eye
%             figure(1); clf; plot(dirs,FRmean,'bo-',dirs,tuning_vonMises(beta,dirs),'r-');
%             pause
            prefDir_fit(n,c) = beta(3);
        end
    end
  
    if plotflag
        figure(n); set(gcf, 'Color', 'w', 'Position', [1 185 370*nCh(n) 770], 'PaperPositionMode', 'auto');
        for c = 1:nCh(n)
            % simple tuning curve
            subplot(2,nCh(n),c);
            dirPlot = dir; dirsPlot = dirs;
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
    
    toc
end

%% next:

% 1) scatter-plot pairwise spike counts/rates across trials of a particular
% stimulus (coh + dir); just do a few pairs to make sure it looks reasonable)
% 2) compute Rsc (noise correlation) of the above, using corr or corrcoef
% 3) plot Rsc as a function of tuning similarity (difference in pref dir)
% later:
% 4) assemble Bondy matrix: binned/smoothed average Rsc as a function of pref dir
% of each member of the pair (will be very sparse until dataset is larger) 





