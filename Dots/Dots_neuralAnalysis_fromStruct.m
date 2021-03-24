% NEW: load/plot data from Miguel struct

clear
close all

% load PrelimDataStruct
load DataStruct_10-26-2020_12-18-2020

plotflag = 1;

extRaster = [50 200]; % extend raster this much before dots onset (1) and after dots offset (2)
latency = [100 150]; % cutoffs for spike rate [vs dots onset (1) and offset (2)]; latency(2) must be < extRaster(2)!
convKernel = fspecial('average', [1 40]); % N ms wide boxcar (acausal, centered)

% some manual plotting lims
psth_xlim1 = 400;
psth_xlim2 = [-300 100];



%% session by session

raster_map = cell(size(dataCell));
spRate_map = cell(size(dataCell));
raster_dotsOn = cell(size(dataCell));
raster_RT = cell(size(dataCell));
spRate_exp = cell(size(dataCell));
spCount_exp = cell(size(dataCell));

for n = 1:length(dataCell)
    n
    %*****************************
    % TUNING (map)
    %*****************************
        
    mStart = dataCell{n}.Mapping.openEvents.motionStart';
        nTrials = length(mStart);
    mEnd = dataCell{n}.Mapping.openEvents.motionEnd';
    dur = round((mEnd-mStart)*1000); % some variability (due to dropped frames???), but okay for now
    dir = dataCell{n}.Mapping.openEvents.direction;
        dirs = unique(dir);
    spd = dataCell{n}.Mapping.openEvents.speed';
        spds = unique(spd);

    apX = nan(size(spd));
    apY = nan(size(spd));
    
    % initialize
    raster_map{n} = nan(length(dataCell{n}.Exp.Channel),nTrials,round(max(dur))+sum(extRaster));
    spRate_map{n} = nan(length(dataCell{n}.Exp.Channel),nTrials); % spike rate on each trial
    
    for t = 1:nTrials 
        apX(t) = mode(dataCell{n}.Mapping.openEvents.apertureX{t}); % or mean/median?
        apY(t) = mode(dataCell{n}.Mapping.openEvents.apertureY{t});
        for c = 1:length(dataCell{n}.Mapping.Channel)
            startTrial = mStart(t)-extRaster(1)/1000;
            endTrial = mEnd(t)+extRaster(2)/1000;
            trialLen = endTrial - startTrial; % in seconds
            raster_map{n}(c,t,1:round(trialLen*1000)) = 0; % init bins to zero
            spkInd = dataCell{n}.Mapping.SpikeTimes{c} >= startTrial & dataCell{n}.Mapping.SpikeTimes{c} <= endTrial;
            spikeTimes = dataCell{n}.Mapping.SpikeTimes{c}(spkInd) - startTrial;
            if ~isempty(spikeTimes)
                spkInd2 = round(spikeTimes*1000);
                if spkInd2(1)==0; spkInd2(1) = 1; end % can't have a zero because used as indices
                raster_map{n}(c,t,spkInd2) = 1;
            end
            a = latency(1)+extRaster(1);
%             b = round(trialLen*1000)-(extRaster(2)-latency(2));
            b = a + dur(t) - latency(1) + latency(2); % same as above but more intuitive
            spRate_map{n}(c,t) = sum(raster_map{n}(c,t,a:b)) / ((b-a+1)/1000);
        end
    end
    pos = [apX apY];
    poss = unique(pos,'rows');

    % fit the tuning curves to get pref dir
    if ~isnan(dataCell{n}.Mapping.SpikeTimes{c})
        for c = 1:length(dataCell{n}.Mapping.Channel)
            FRmean = calc_mean(spRate_map{n}(c,:)',dir,dirs);
            peak = dirs(FRmean==max(FRmean));
            % vonMises params: b = [A, k, theta, offset]
            guess = [max(FRmean)-min(FRmean), 1.5, peak(1), min(FRmean)];
            tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
            tuning_vonMises_err = @(b,dir,FR) sum((tuning_vonMises(b,dir)-FR).^2);
            [beta, fval] = fminsearch(@(x) tuning_vonMises_err(x,dir',spRate_map{n}(c,:)), guess);
    %        % temp, to check fits by eye: (they're all good!)
    %         err = fval/length(dir)
    %         figure(1); clf; plot(dirs,FRmean,'bo-',dirs,tuning_vonMises(beta,dirs),'r-');
    %         pause
            dataCell{1}.Mapping.Selectivity_fit(c) = beta(3);
        end
    end
    
    if plotflag
    figure(n); set(gcf, 'Color', 'w', 'Position', [200 500 370*length(dataCell{n}.Exp.Channel) 770], 'PaperPositionMode', 'auto');
    for c = 1:length(dataCell{n}.Mapping.Channel)
        % simple tuning curve
        subplot(2,length(dataCell{n}.Mapping.Channel),c);
        [FRmean, FRse, ~] = calc_mean(spRate_map{n}(c,:)',dir,dirs);
%         FRmean(end+1) = FRmean(1); FRse(end+1) = FRse(1); dirs3 = [dirs 360]; % kluge for standard set of dirs, be careful
        errorbar(dirs,FRmean,FRse,'ko-','MarkerFaceColor','k');
        set(gca,'xtick',dirs,'xticklabel',{'0','','90','','180','','270',''},'tickdir','out')
        xlabel('direction (deg)'); ylabel('firing rate (sp/s)');
        title(sprintf('sess %d, chan %d',n,dataCell{n}.Mapping.Channel(c)));
        changeAxesFontSize(gca, 22, 22);
        
        % could repeat for speed, pos, size.. but only really needed during the expt
        
        % PSTH for pref/null, aligned to dots onset
        subplot(2,length(dataCell{n}.Mapping.Channel),c+length(dataCell{n}.Mapping.Channel));
        pref = dirs(FRmean==max(FRmean)); null = dirs(FRmean==min(FRmean));
        psth = calc_mean(squeeze(raster_map{n}(c,:,:)), dir, [pref;null])*1e3;
        psth = smoothRaster(psth, convKernel);
%         tAxis = -extRaster(1) : mode(dur)+latency(2);
        tAxis = -extRaster(1) : size(squeeze(raster_map{n}(c,:,:)),2) - extRaster(1) - 1;
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
%         title(sprintf('sess %d, chan %d',n,dataCell{n}.Mapping.Channel(c)));
        changeAxesFontSize(gca, 22, 22);
    end
    end
    
    %*****************************
    % TASK (exp)
    %*****************************
    
    mStart = dataCell{n}.Exp.openEvents.motionStart';
        nTrials = length(mStart);
    mEnd = dataCell{n}.Exp.openEvents.motionEnd';
    endAcq = dataCell{n}.Exp.openEvents.endAcquire';
    RT = round((mEnd-mStart)*1000);
%     RT = round((endAcq-mStart)*1000); % if used, this would need to subtract the 'hold' requirement (but I don't think it should be used)
    dur=RT;
    dir = dataCell{n}.Exp.openEvents.direction;
        dirs = unique(dir);

    % initialize
    raster_dotsOn{n} = nan(length(dataCell{n}.Exp.Channel),nTrials,round(max(dur))+sum(extRaster));
    raster_RT{n} =     nan(length(dataCell{n}.Exp.Channel),nTrials,round(max(dur))+sum(extRaster));
    spRate_exp{n} =    nan(length(dataCell{n}.Exp.Channel),nTrials);
    spCount_exp{n} =   nan(length(dataCell{n}.Exp.Channel),nTrials);

    for t = 1:nTrials    
        for c = 1:length(dataCell{n}.Exp.Channel)
            startTrial = mStart(t)-extRaster(1)/1000;
            endTrial = mEnd(t)+extRaster(2)/1000;
            trialLen = endTrial - startTrial; % in seconds
            
            % raster, align dotsOn
            raster_dotsOn{n}(c,t,1:round(trialLen*1000)) = 0; % init bins to zero but keep nans after trial end
            spkInd = dataCell{n}.Exp.SpikeTimes{c} >= startTrial & dataCell{n}.Exp.SpikeTimes{c} <= endTrial;
            spikeTimes = dataCell{n}.Exp.SpikeTimes{c}(spkInd) - startTrial;
            if ~isempty(spikeTimes)
                spkInd2 = round(spikeTimes*1000);
                if spkInd2(1)==0; spkInd2(1) = 1; end % can't have a zero because used as indices
                raster_dotsOn{n}(c,t,spkInd2) = 1;
            end
            
            % raster, align RT (dotsOff)
            if round(trialLen*1000) == size(raster_RT{n},3)
                raster_RT{n}(c,t,:) = 0;
            else
                raster_RT{n}(c,t,end-round(trialLen*1000):end) = 0; % init bins to zero but keep nans before dots on
            end
            spikeTimes = -(endTrial - dataCell{n}.Exp.SpikeTimes{c}(spkInd)); % minus because counting backwards from endTrial
            if ~isempty(spikeTimes)
                spkInd2 = round(spikeTimes*1000) + size(squeeze(raster_RT{n}(c,:,:)),2); % plus here because the minus was already applied
                if spkInd2(1)==0; spkInd2(1) = 1; end % can't have a zero because used as indices
                if spkInd2(end)>size(raster_RT{n},3); spkInd2(end) = size(raster_RT{n},3); end
                raster_RT{n}(c,t,spkInd2) = 1;
            end
            
            % spike rate
            a = latency(1)+extRaster(1);
%             b = round(trialLen*1000)-(extRaster(2)-latency(2));
            b = a + dur(t) - latency(1) + latency(2); % same as above but more intuitive
            spRate_exp{n}(c,t) = sum(raster_dotsOn{n}(c,t,a:b)) / ((b-a+1)/1000);
            spCount_exp{n}(c,t) = sum(raster_dotsOn{n}(c,t,a:b));
        end
    end

    if plotflag
    % choice+conf-conditioned PSTHs
    J = abs(dataCell{n}.Exp.openEvents.coherence)==0; % weak motion only, or else there's a confound! can also use residuals...

        % aligned dots-on
    figure(n*10); set(gcf, 'Color', 'w', 'Position', [200 300 370*length(dataCell{n}.Exp.Channel) 430], 'PaperPositionMode', 'auto');
    for c = 1:length(dataCell{n}.Exp.Channel)
        subplot(1,length(dataCell{n}.Mapping.Channel),c);
        clear psth;
        I = J & dataCell{n}.Exp.openEvents.choice==1 & dataCell{n}.Exp.openEvents.pdw==1; % left-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_dotsOn{n}(c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==1 & dataCell{n}.Exp.openEvents.pdw==0; % left-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_dotsOn{n}(c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==2 & dataCell{n}.Exp.openEvents.pdw==1; % right-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_dotsOn{n}(c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==2 & dataCell{n}.Exp.openEvents.pdw==0; % right-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_dotsOn{n}(c,I,:)))*1e3, convKernel);
        
        tAxis = -extRaster(1) : min([1000 size(raster_dotsOn{n},3)-sum(extRaster)]);
        tAxis = tAxis + length(convKernel)/2+1; % offset for causal
        clr = {'g-', 'g--', 'b-', 'b--'};
        YLim = [min(min(psth))*0.9 max(max(psth))*1.1];
        XLim = [tAxis(1) psth_xlim1];
        hold on;
        for p = 1:size(psth,1)
            plot(tAxis, psth(p,1:length(tAxis)), clr{p},'LineWidth',2);
        end
        set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
        xlabel('Time from dots on (ms)');
        ylabel('Firing rate (spikes/s)');
        changeAxesFontSize(gca, 22, 22);
        if c==length(dataCell{n}.Exp.Channel)
            legend('left high','left low','right high','right low','Location','Northeast'); legend('boxoff');
        end
    end

        % aligned RT
    figure(n*100); set(gcf, 'Color', 'w', 'Position', [200 100 370*length(dataCell{n}.Exp.Channel) 430], 'PaperPositionMode', 'auto');
    for c = 1:length(dataCell{n}.Exp.Channel)
        subplot(1,length(dataCell{n}.Mapping.Channel),c);
        clear psth
        I = J & dataCell{n}.Exp.openEvents.choice==1 & dataCell{n}.Exp.openEvents.pdw==1; % left-high
        psth(1,:) = smoothRaster(nanmean(squeeze(raster_RT{n}(c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==1 & dataCell{n}.Exp.openEvents.pdw==0; % left-low
        psth(2,:) = smoothRaster(nanmean(squeeze(raster_RT{n}(c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==2 & dataCell{n}.Exp.openEvents.pdw==1; % right-high
        psth(3,:) = smoothRaster(nanmean(squeeze(raster_RT{n}(c,I,:)))*1e3, convKernel);
        I = J & dataCell{n}.Exp.openEvents.choice==2 & dataCell{n}.Exp.openEvents.pdw==0; % right-low
        psth(4,:) = smoothRaster(nanmean(squeeze(raster_RT{n}(c,I,:)))*1e3, convKernel);
        
        tAxis = -(size(raster_RT{n},3)-extRaster(2)-1) : extRaster(2);
        tAxis = tAxis + length(convKernel)/2+1; % offset for causal
        XLim = psth_xlim2;
        YLim = [min(min(psth))*0.9 max(max(psth))*1.1];
        hold on;
        for p = 1:size(psth,1)
            plot(tAxis, psth(p,:), clr{p},'LineWidth',2);
        end
        set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
        xlabel('Time from saccade (ms)');
        changeAxesFontSize(gca, 22, 22);
    end
    
    end
    
    % next: aggregate data across cells
    
    

    
end












