function dots3DMP_plots_fit(Data,fitInterp,mods,cohs,deltas,hdgs,conftask,RTtask,fitgauss) %#ok<INUSL>

% somewhat unwieldy amalgam of parseData and dots3DMP_plots, to show model 
% fits versus data in dots3DMP experiment

% calls parseData separately on fit and data, then plots the curves for
% the former and just the data points for the latter

if nargin<3, conftask=1; end
if nargin<4, RTtask = 0; end

%% data

data = Data;
mods = unique(data.modality);
cohs = unique(data.coherence);
deltas = unique(data.delta);
hdgs = unique(data.heading);

dots3DMP_parseData

% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR omit the zero delta, for some plots
%     D = find(deltas~=0)';

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c); %#ok<*NODEF>
    for m = 1:length(mods)     % m c d h
        h(m) = errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), [clr{c}{m}]); hold on; %#ok<*IDISVAR>
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));
    for m = 1:length(mods)
        h(m) = errorbar(hdgs, squeeze(confMean(m,c,D,:)), squeeze(confSE(m,c,D,:)), [clr{c}{m}]);
        ylim([0 1]); hold on;
    end
    xlabel('heading angle (deg)'); 
    if conftask==1, ylabel('saccadic endpoint (''confidence'', %)');
    elseif conftask==2, ylabel('proportion high bet');
    end
    
    if RTtask
        subplot(3,length(cohs),c+length(cohs)*2);
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:)), squeeze(RTse(m,c,D,:)), [clr{c}{m}]); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end

% now separate by delta
if length(deltas)>1

clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
figure(108);
set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2+RTtask,length(cohs),c);
    for d = 1:length(deltas)     % m c d h
        h(d) = errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), [clr{c}{d}]); hold on
        L{d} = sprintf('\x0394=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));
    for d = 1:length(deltas)
        h(d) = errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), [clr{c}{d}]);
        ylim([0 1]); hold on;
    end
    xlabel('heading angle (deg)'); 
    if conftask==1, ylabel('saccadic endpoint (''confidence'', %)');
    elseif conftask==2, ylabel('proportion high bet');
    end  
    
    if RTtask
        subplot(3,length(cohs),c+length(cohs)*2);
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(RTmean(3,c,d,:)), squeeze(RTse(3,c,d,:)), [clr{c}{d}]); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end

end


%% fit

data = fitInterp;

mods = unique(data.modality);
cohs = unique(data.coherence);
deltas = unique(data.delta);
hdgs = unique(data.heading);

dots3DMP_parseData

if fitgauss
    
    dots3DMP_fit_cgauss
    
    exfig=0;

    D = 1:length(deltas);

% % % %     clr{1} = {[1 0 0], [0 0 1]}; % red, blue
% % % %     clr{1} = {[1 0 1], [0 1 1]}; % magenta, cyan (TMS)
% % %     clr{1} = {[0 0 1], [1 0 0]}; % blue (baseline), red (TMS)
% % %     clr{2} = clr{1};
% % %     clr{3} = clr{1};

    for d = D
        % choice
        clear L h; figure(300+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
        lind = 1;
        for c = 1:length(cohs)
            beta = [muPMF(3,c,d) sigmaPMF(3,c,d)];
            h(lind) = plot(xVals, cgauss(beta,xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
%             errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
            L{lind} = sprintf('Coh=%0.1f',cohs(c)); lind = lind+1;
            xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
            changeAxesFontSize(gca,20,20); set(gca,'box','off')
            legend(h,L,'location','northwest'); legend('boxoff');
        end
        if exfig; export_fig(['d=' num2str(deltas(d)) '_pmf'],'-eps'); end

        % conf
        figure(400+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
        for c = 1:length(cohs)
            beta = [amplConf(3,c,d) muConf(3,c,d) sigmaConf(3,c,d) baselineConf(3,c,d)];
            plot(xVals, flippedGauss(beta,xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
%             errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
            xlabel('Heading angle (deg)'); ylabel('Confidence'); ylim([0 1]);
            changeAxesFontSize(gca,20,20); set(gca,'box','off')
        end
        if exfig; export_fig(['d=' num2str(deltas(d)) '_conf'],'-eps'); end

        % RT
        if ~isnan(RTmean(1,1,2))
            figure(500+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
            for c = 1:length(cohs)
                beta = [amplRT(3,c,d) muRT(3,c,d) sigmaRT(3,c,d) baselineRT(3,c,d)];
                plot(xVals, gauss(beta,xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
%                 errorbar(hdgs, squeeze(RTmean(3,c,d,:)), squeeze(RTse(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
                xlabel('Heading angle (deg)'); ylabel('Reaction time (s)');
                changeAxesFontSize(gca,20,20); set(gca,'box','off')
            end
            if exfig; export_fig(['d=' num2str(deltas(d)) '_rt'],'-eps'); end
        end    
    end

    
else
         %ves %vis %comb
clr{1} = {'k-','m-','c-'};
clr{2} = {'k-','r-','b-'};
clr{3} = {'k-','y-','g-'};
figure(101+D);
for c = 1:length(cohs)
    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c);
    for m = 1:length(mods)     % m c d h
        h(m) = plot(hdgs, squeeze(pRight(m,c,D,:)), [clr{c}{m}]); hold on;
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));
    for m = 1:length(mods)
        h(m) = plot(hdgs, squeeze(confMean(m,c,D,:)), [clr{c}{m}]);
        ylim([0 1]); hold on;
    end
    xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');
    
    if ~isnan(RTmean(1,1,2))
        subplot(3,length(cohs),c+length(cohs)*2);
        for m = 1:length(mods)
            h(m) = plot(hdgs, squeeze(RTmean(m,c,D,:)), [clr{c}{m}]); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end


% now separate by delta

if length(deltas)>1

clr{1} = {'b-','c-','g-'};
clr{2} = {'b-','c-','g-'};
clr{3} = {'b-','c-','g-'};

clear L;
figure(108);
for c = 1:length(cohs)
    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c);
    for d = 1:length(deltas)     % m c d h
        h(d) = plot(hdgs, squeeze(pRight(3,c,d,:)), [clr{c}{d}]); hold on
        L{d} = sprintf('\x0394=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));
    for d = 1:length(deltas)
        h(d) = plot(hdgs, squeeze(confMean(3,c,d,:)), [clr{c}{d}]);
        ylim([0 1]); hold on;
    end
    xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');
    
    if ~isnan(RTmean(1,1,2))
        subplot(3,length(cohs),c+length(cohs)*2);
        for d = 1:length(deltas)
            h(d) = plot(hdgs, squeeze(RTmean(3,c,d,:)), [clr{c}{d}]); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end

end

end



