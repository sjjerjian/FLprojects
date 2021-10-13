function dots3DMP_plots_fit_byCoh(rawData,fitInterp,conftask,RTtask,fitgauss)

% somewhat unwieldy amalgam of parseData and dots3DMP_plots, to show model 
% fits versus data in dots3DMP experiment

% calls parseData separately on fit and data, then plots the curves for
% the former and just the data points for the latter
% if fitgauss == 1, fit curves will be gaussian fits, otherwise they will
% just be interpolated between headings

if nargin<3, conftask=1; end
if nargin<4, RTtask = 0; end
if nargin<5, fitgauss = 0; end

%% actual data

mods   = unique(rawData.modality);
cohs   = unique(rawData.coherence);
deltas = unique(rawData.delta);
hdgs   = unique(rawData.heading);

parsedData = dots3DMP_parseData(rawData,mods,cohs,deltas,hdgs,conftask,RTtask); 

spRows = 1 + double(conftask>0) + double(RTtask);

if conftask==1, confYL = 'Sacc EP';
elseif conftask==2, confYL = 'P(High Bet)';
end

D = 2;
% first, for all trials irrespective of delta
% D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR omit the zero delta, for some plots
%     D = find(deltas~=0)';

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c); %#ok<*NODEF>
    for m = 1:length(mods)     % m c d h
        h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), [clr{c}{m}]); hold on; %#ok<*IDISVAR>
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('P(Right)');

    if conftask>0
        subplot(spRows,length(cohs),c+length(cohs));
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), [clr{c}{m}]);
            ylim([0 1]); hold on;
        end
        xlabel('heading angle (deg)');
        ylabel(confYL);
    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0)));
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), [clr{c}{m}]); hold on;
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
    subplot(spRows,length(cohs),c);
    for d = 1:length(deltas)     % m c d h
        h(d) = errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), [clr{c}{d}]); hold on
        L{d} = sprintf('\x0394=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');
    
    if conftask>0
        subplot(spRows,length(cohs),c+length(cohs));
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), [clr{c}{d}]);
            ylim([0 1]); hold on;
        end
        xlabel('heading angle (deg)');
        ylabel(confYL);
    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0)));
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), [clr{c}{d}]); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end

end


%% fit

% parsedData = dots3DMP_parseData(fitInterp,mods,cohs,deltas,hdgs,conftask,RTtask); 

mods   = unique(fitInterp.modality);
cohs   = unique(fitInterp.coherence);
deltas = unique(fitInterp.delta);
hdgs   = unique(fitInterp.heading);

if fitgauss
    
%     gfit = dots3DMP_fit_cgauss(fitInterp,mods,cohs,deltas,conftask,RTtask);
    
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
            beta = [gfit.choice.mu(3,c,d) gfit.choice.sigma(3,c,d)];
            h(lind) = plot(parsedData.xVals, cgauss(beta,parsedData.xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
%             errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
            L{lind} = sprintf('Coh=%0.1f',cohs(c)); lind = lind+1;
            xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
            changeAxesFontSize(gca,20,20); set(gca,'box','off')
            legend(h,L,'location','northwest'); legend('boxoff');
        end
        if exfig; export_fig(['d=' num2str(deltas(d)) '_pmf'],'-eps'); end

        % conf
        if conftask > 0
            figure(400+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
            for c = 1:length(cohs)
                beta = [gfit.conf.ampl(3,c,d) gfit.conf.mu(3,c,d) gfit.conf.sigma(3,c,d) gfit.conf.bsln(3,c,d)];
                plot(parsedData.xVals, flippedGauss(beta,parsedData.xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
                %             errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
                xlabel('Heading angle (deg)');
                ylabel(confYL);
                ylim([0 1]);
                changeAxesFontSize(gca,20,20); set(gca,'box','off')
            end
            if exfig; export_fig(['d=' num2str(deltas(d)) '_conf'],'-eps'); end
        end
        
        % RT
        if RTtask
            figure(500+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
            for c = 1:length(cohs)
                beta = [gfit.RT.ampl(3,c,d) gfit.RT.mu(3,c,d) gfit.RT.sigma(3,c,d) gfit.RT.bsln(3,c,d)];
                plot(parsedData.xVals, gauss(beta,parsedData.xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
                %                 errorbar(hdgs, squeeze(RTmean(3,c,d,:)), squeeze(RTse(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
                xlabel('Heading angle (deg)'); ylabel('RT (s)');
                changeAxesFontSize(gca,20,20); set(gca,'box','off')
            end
            if exfig; export_fig(['d=' num2str(deltas(d)) '_rt'],'-eps'); end
        end
    end

    
else
    
D = 2;
         %ves %vis %comb
clr{1} = {'k-','m-','c-'};
clr{2} = {'k-','r-','b-'};
clr{3} = {'k-','y-','g-'};
figure(101+D);
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c);
    for m = 1:length(mods)     % m c d h
        h(m) = plot(hdgs, squeeze(fitInterp.pRight(m,c,D,:)), [clr{c}{m}]); hold on;
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('P(Right)');

    if conftask>0
    subplot(spRows,length(cohs),c+length(cohs));
    for m = 1:length(mods)
        h(m) = plot(hdgs, squeeze(fitInterp.confMean(m,c,D,:)), [clr{c}{m}]);
        ylim([0 1]); hold on;
    end
    xlabel('heading angle (deg)'); 
    ylabel(confYL);
    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*2);
        for m = 1:length(mods)
            h(m) = plot(hdgs, squeeze(fitInterp.RTmean(m,c,D,:)), [clr{c}{m}]); hold on;
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
    subplot(spRows,length(cohs),c);
    for d = 1:length(deltas)     % m c d h
        h(d) = plot(hdgs, squeeze(fitInterp.pRight(3,c,d,:)), [clr{c}{d}]); hold on
        L{d} = sprintf('\x0394=%d',deltas(d));
    end
    if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    ylim([0 1]);
    
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('P(Right)');
    
    if conftask>0
        subplot(spRows,length(cohs),c+length(cohs));
        for d = 1:length(deltas)
            h(d) = plot(hdgs, squeeze(fitInterp.confMean(3,c,d,:)), [clr{c}{d}]);
            hold on;
        end
        xlabel('heading angle (deg)'); ylabel(confYL);
        ylim([0 1]); 
    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0)));
        for d = 1:length(deltas)
            h(d) = plot(hdgs, squeeze(fitInterp.RTmean(3,c,d,:)), [clr{c}{d}]); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end

end

end



