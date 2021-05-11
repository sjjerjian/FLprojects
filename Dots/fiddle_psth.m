% fiddle_pseudopop.m

clear
close all
load('/Users/chris/Documents/MATLAB/Projects/offlineTools/Dots/PseudoPopulationOne.mat')
data = data_;
clear data_


offset = 50; % add offset (ms) to the front and back end of each raster
convKernel = fspecial('average', [1 40]); % N ms wide boxcar (acausal, centered)

% manual axis lims for plotting
psth_xmax = 400;
psth_xmax2 = 200;
YMax(11) = 250;
YMax(21) = 100;
YMax(31) = 350;
YMax(32) = 400;
YMax(33) = 80;

% some cleanup
for s = 1:length(data)
    data{s}.spikeTime(data{s}.spikeTime==0) = NaN;
    data{s} = rmfield(data{s},'timeDirection');
end

% keep only 12, 21, and 4x
data{1}.Channel([1 3]) = [];
data{1}.spikeTime([1 3],:) = [];
data{2}.Channel(2) = [];
data{2}.spikeTime(2,:) = [];

% pool 4 + 5
addTime = data{4}.endAcquire(end)+10; % end of block 4 plus 10 sec, add this to block 5's times
data{5}.motionStart = data{5}.motionStart + addTime;
data{5}.motionEnd = data{5}.motionEnd + addTime;
data{5}.FixPoint = data{5}.FixPoint + addTime;
data{5}.FPAcquire = data{5}.FPAcquire + addTime;
data{5}.endAcquire = data{5}.endAcquire + addTime;
data{5}.trialStart = data{5}.trialStart + addTime;
data{5}.spikeTime = data{5}.spikeTime + addTime;
data{5}.trialNum = data{5}.trialNum + data{4}.trialNum(end); % same idea w trial num
fn = fieldnames(data{4});
for F = 1:length(fn)
    if ~strcmp(fn{F},'Channel') && ~strcmp(fn{F},'Date')
        eval(['data{4}.' fn{F} ' = [data{4}.' fn{F} ' data{5}.' fn{F} '];']);
    end
end

% get rid of 5 and 3, 4 becomes 3
data = data([1 2 4]);

for s = 2 %:length(data) % just do one at a time for demo purposes
   
%     % group cohs for simpler plots
%     data{s}.coherence(data{s}.coherence>0.25) = 2;
%     data{s}.coherence(data{s}.coherence<0.25 & data{s}.coherence>0) = 1;
%     % or:
% %     data{s}.coherence(data{s}.coherence>0) = 1;

    % convert to signed coh
    if any(~ismember(data{s}.direction,[0 180]))
        keyboard
    end
    data{s}.coherence(data{s}.direction==180) = -data{s}.coherence(data{s}.direction==180);

    cohs = unique(data{s}.coherence);
%     cohs(cohs==0) = []; % leave out 0 coh in some plots?

    % some useful session vars
    ntrials = length(data{s}.motionStart);
%     RT = data{s}.endAcquire-data{s}.motionStart-0.08; % 0.08 because 80 ms from target acquired to endAcquire [double check that]
    RT = data{s}.motionEnd-data{s}.motionStart; % motionEnd is better!

    raster_dotsOn = cell(length(data{s}.Channel),1);
    raster_RT = cell(length(data{s}.Channel),1);
    for n = 1:length(data{s}.Channel)
        
        % initialize
        raster_dotsOn{n} = nan(ntrials,round(max(RT)*1000)+offset*2); % add offset to the front and back end of each raster, so total of offset*2
        raster_RT{n} = nan(ntrials,round(max(RT)*1000)+offset*2);
        
        for t = 1:ntrials
            startTrial = data{s}.motionStart(t)-offset/1000;
            endTrial = data{s}.motionEnd(t)+offset/1000;
            trialLen = endTrial - startTrial;
            
            % align dotsOn
            raster_dotsOn{n}(t,1:round(trialLen*1000)) = 0; % init bins to zero but keep nans after trial end
            spkInd = data{s}.spikeTime(n,:) >= startTrial & data{s}.spikeTime(n,:) <= endTrial;
            spikeTimes = data{s}.spikeTime(n,spkInd) - startTrial;
            if ~isempty(spikeTimes)
                spkInd2 = round(spikeTimes*1000);
                if spkInd2(1)==0; spkInd2(1) = 1; end
                raster_dotsOn{n}(t,spkInd2) = 1;
            end
            
            % align RT (dotsOff)
            if round(trialLen*1000) == size(raster_RT{n},2)
                raster_RT{n}(t,:) = 0;
            else
                raster_RT{n}(t,end-round(trialLen*1000):end) = 0; % init bins to zero but keep nans before dots on
            end
            spikeTimes = -(endTrial - data{s}.spikeTime(n,spkInd)); % minus because counting backwards from endTrial
            if ~isempty(spikeTimes)
                spkInd2 = round(spikeTimes*1000) + size(raster_RT{n},2); % plus here because the minus was already applied
                if spkInd2(1)==0; spkInd2(1) = 1; end
                if spkInd2(end)>size(raster_RT{n},2); spkInd2(end) = size(raster_RT{n},2); end
                raster_RT{n}(t,spkInd2) = 1;
            end
        end
        
        % coh-conditioned PSTH, aligned to dots onset and RT
        figure(s*10+n); set(gcf, 'Color', 'w', 'Position', [100+20*(s*10+n) 800 500 500], 'PaperPositionMode', 'auto');
            % aligned dots-on
        psth = calc_mean(raster_dotsOn{n},data{s}.coherence', cohs')*1e3;
        psth = smoothRaster(psth, convKernel);
        tAxis = -offset : min([1000 size(raster_dotsOn{n},2)]);
        tAxis = tAxis + length(convKernel)/2+1; % offset for causal
        clr = cool(length(cohs));
        YLim = [min(min(psth))*0.9 max(max(psth))*1.1];
        XLim = [tAxis(1) psth_xmax];
        hold on;
        for c = 1:length(cohs)
            h(c) = plot(tAxis, psth(c,1:length(tAxis)),'Color', clr(c,:),'LineWidth',2);
        end
        set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
        xlabel('Time from dots on (ms)');
        ylabel('Firing rate (spikes/s)');
        ylim([0 YMax(s*10+n)]);
        changeAxesFontSize(gca, 22, 22);
        legend([h(1) h(4) h(8) h(11)],'strong left motion','weak left motion','weak right motion','strong right motion','Location','Northeast'); legend('boxoff');
              % ^ hard coded for 11 cohs, can try to make it more general?
%         keyboard % to move legend before exporting fig
%         export_fig([num2str(s*10+n) '_cohs_A'], '-eps');
        
            % aligned RT
        figure(s*10+n+1000); set(gcf, 'Color', 'w', 'Position', [600+20*(s*10+n) 800 335 500], 'PaperPositionMode', 'auto');
        psth = calc_mean(raster_RT{n},data{s}.coherence', cohs')*1e3;
        psth = smoothRaster(psth, convKernel);
        tAxis = -(size(raster_RT{n},2)-offset-1) : offset;
        tAxis = tAxis + length(convKernel)/2+1; % offset for causal
        XLim = [-psth_xmax2 tAxis(end)];
        hold on;
        for c = 1:length(cohs)
            h(c) = plot(tAxis, psth(c,:),'Color', clr(c,:),'LineWidth',2);
        end
        set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
        xlabel('Time from saccade (ms)');
        ylim([0 YMax(s*10+n)]);
        changeAxesFontSize(gca, 22, 22);
%         export_fig([num2str(s*10+n) '_cohs_B'], '-eps');
        
        
        %%%%%%%%%
        
        % conf-conditioned PSTHs
        figure(s*100+n); set(gcf, 'Color', 'w', 'Position', [100+20*(s*10+n) 200 500 500], 'PaperPositionMode', 'auto');
        
            % aligned dots-on
        clear psth;
        J = abs(data{s}.coherence)<=0.032; % weak motion only, or else there's a confound; can also use residuals!
        I = J & data{s}.choice==1 & data{s}.pdw==1; % left-high
        psth(1,:) = smoothRaster(nanmean(raster_dotsOn{n}(I,:))*1e3, convKernel);
        I = J & data{s}.choice==1 & data{s}.pdw==0; % left-low
        psth(2,:) = smoothRaster(nanmean(raster_dotsOn{n}(I,:))*1e3, convKernel);
        I = J & data{s}.choice==2 & data{s}.pdw==1; % right-high
        psth(3,:) = smoothRaster(nanmean(raster_dotsOn{n}(I,:))*1e3, convKernel);
        I = J & data{s}.choice==2 & data{s}.pdw==0; % right-low
        psth(4,:) = smoothRaster(nanmean(raster_dotsOn{n}(I,:))*1e3, convKernel);
        
        tAxis = -offset : min([1000 size(raster_dotsOn{n},2)]);
        tAxis = tAxis + length(convKernel)/2+1; % offset for causal
        clr = {'g-', 'g--', 'b-', 'b--'};
        YLim = [min(min(psth))*0.9 max(max(psth))*1.1];
        XLim = [tAxis(1) psth_xmax];
        hold on;
        for c = 1:size(psth,1)
            h(c) = plot(tAxis, psth(c,1:length(tAxis)), clr{c},'LineWidth',2);
        end
        set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out'); box off;
        xlabel('Time from dots on (ms)');
        ylabel('Firing rate (spikes/s)');
        ylim([0 YMax(s*10+n)]);
        changeAxesFontSize(gca, 22, 22);
        legend('left high','left low','right high','right low','Location','Northeast'); legend('boxoff');
%         export_fig([num2str(s*10+n) '_confP_A'], '-eps');

            % aligned RT
        figure(s*100+n+1000); set(gcf, 'Color', 'w', 'Position', [600+20*(s*10+n) 200 335 500], 'PaperPositionMode', 'auto');
        clear psth
        I = J & data{s}.choice==1 & data{s}.pdw==1; % left-high
        psth(1,:) = smoothRaster(nanmean(raster_RT{n}(I,:))*1e3, convKernel);
        I = J & data{s}.choice==1 & data{s}.pdw==0; % left-low
        psth(2,:) = smoothRaster(nanmean(raster_RT{n}(I,:))*1e3, convKernel);
        I = J & data{s}.choice==2 & data{s}.pdw==1; % right-high
        psth(3,:) = smoothRaster(nanmean(raster_RT{n}(I,:))*1e3, convKernel);
        I = J & data{s}.choice==2 & data{s}.pdw==0; % right-low
        psth(4,:) = smoothRaster(nanmean(raster_RT{n}(I,:))*1e3, convKernel);
        
        tAxis = -(size(raster_RT{n},2)-offset-1) : offset;
        tAxis = tAxis + length(convKernel)/2+1; % offset for causal
        XLim = [-psth_xmax2 tAxis(end)];
        hold on;
        for c = 1:size(psth,1)
            h(c) = plot(tAxis, psth(c,:), clr{c},'LineWidth',2);
        end
        set(gca, 'XLim', XLim, 'YLim', YLim, 'TickDir','out', 'YAxisLocation', 'Right'); box off;
        xlabel('Time from saccade (ms)');
        ylim([0 YMax(s*10+n)]);
        changeAxesFontSize(gca, 22, 22);
%         export_fig([num2str(s*10+n) '_confP_B'], '-eps');

    end
end

         
      