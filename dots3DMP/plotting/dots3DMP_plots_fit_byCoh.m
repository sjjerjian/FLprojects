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

fsz = 14; % fontsize

mods   = unique(rawData.modality);
cohs   = unique(rawData.coherence);
deltas = unique(rawData.delta);
hdgs   = unique(rawData.heading);

if all(mods==1), cohs=1; end

xLab = sprintf('heading angle (%s)',char(176));
modlabels = {'Ves','Vis','Comb'};

if conftask==1 
    yLab = 'SaccEP'; confYlims = [0.2 0.9]; 
    if all(mods==1), RTylims = [0.9 1.5]; 
    else,            RTylims = [0.9 1.9]; 
    end
    xt = -10:5:10;
     if length(mods)>1, cohlabs = {'Low Coh','High Coh'}; end
elseif conftask==2
    yLab = 'P(High Bet)'; confYlims = [0.4 1.0]; 
    if all(mods==1), RTylims = [0.5 0.72]; 
    else,            RTylims = [0.5 0.9]; 
    end
    xt = -12:6:12;
    if length(mods)>1, cohlabs = {['coh = ' num2str(cohs(1))],['coh = ' num2str(cohs(2))]}; end
end 

%% actual data
% plot as individual points with error bars

useAbsHdg = 0;
parsedData = dots3DMP_parseData(rawData,mods,cohs,deltas,hdgs,conftask,RTtask,useAbsHdg); 

if useAbsHdg, hdgs = unique(abs(hdgs)); end

spRows = 1 + double(conftask>0) + double(RTtask);

D = find(deltas==0);
% first, for all trials irrespective of delta
% D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR omit the zero delta, for some plots
%     D = find(deltas~=0)';

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 1000 230+300*(length(cohs)-1) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c); hold on;
    for m = 1:length(mods)     % m c d h
        h(m) = errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), [clr{c}{m}],'linewidth',1.5); 
    end
    
    if conftask>0
        subplot(spRows,length(cohs),c+length(cohs)); hold on
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), [clr{c}{m}],'linewidth',1.5);
        end

    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); hold on
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), [clr{c}{m}],'linewidth',1.5); 
        end
    end
    
end

% now separate by delta
if length(deltas)>1

clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
figure(108);
set(gcf,'Color',[1 1 1],'Position',[50 20 230+300*(length(cohs)-1) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(spRows,length(cohs),c); hold on
    for d = 1:length(deltas)     % m c d h
        h(d) = errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), [clr{c}{d}],'linewidth',1.5);
%         L{d} = sprintf('\x0394=%d',deltas(d));
    end

    if conftask>0
        subplot(spRows,length(cohs),c+length(cohs)); hold on
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), [clr{c}{d}],'linewidth',1.5);
        end
    end
    
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); hold on
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), [clr{c}{d}],'linewidth',1.5); 
        end
    end
    
end

end


%% fit
% plot curve only, overlaid on actual data points

% parsedData = dots3DMP_parseData(fitInterp,mods,cohs,deltas,hdgs,conftask,RTtask); 

mods   = unique(fitInterp.modality);
cohs   = unique(fitInterp.coherence);
deltas = unique(fitInterp.delta);
hdgs   = unique(fitInterp.heading);
D = find(deltas==0);

if fitgauss
    
%     gfit = dots3DMP_fit_cgauss(fitInterp,mods,cohs,deltas,conftask,RTtask);
    
    exfig=0;
    clear L h; 

% % % %     clr{1} = {[1 0 0], [0 0 1]}; % red, blue
% % % %     clr{1} = {[1 0 1], [0 1 1]}; % magenta, cyan (TMS)
% % %     clr{1} = {[0 0 1], [1 0 0]}; % blue (baseline), red (TMS)
% % %     clr{2} = clr{1};
% % %     clr{3} = clr{1};
    
    figure(108)
    for d = 1:length(deltas)
        % choice
        lind = 1;
        for c = 1:length(cohs)
            beta = [gfit.choice.mu(3,c,d) gfit.choice.sigma(3,c,d)];
            h(lind) = plot(parsedData.xVals, cgauss(beta,parsedData.xVals), '-o', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
%             errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
            L{lind} = sprintf('Coh=%0.1f',cohs(c)); lind = lind+1;
            xlabel('Heading angle (deg)'); ylabel('Proportion rightward choices'); ylim([0 1]);
            changeAxesFontSize(gca,20,20); set(gca,'box','off')
            legend(h,L,'location','northwest'); legend('boxoff');
        end
        if length(mods)>1; title(cohlabs{c}); end
        
        % conf
        if conftask > 0
            figure(400+d); set(gcf,'Color',[1 1 1],'Position',[50 20 360 320],'PaperPositionMode','auto'); clf;
            for c = 1:length(cohs)
                beta = [gfit.conf.ampl(3,c,d) gfit.conf.mu(3,c,d) gfit.conf.sigma(3,c,d) gfit.conf.bsln(3,c,d)];
                plot(parsedData.xVals, flippedGauss(beta,parsedData.xVals), '-', 'Color', clr{d}{c}, 'Linewidth', 3); hold on;
                %             errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), 'o', 'Color', clr{d}{c}, 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2);
                xlabel(xLab);
                ylabel(yLab);
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
    
    %ves %vis %comb
    clr{1} = {'k-','m-','c-'};
    clr{2} = {'k-','r-','b-'};
    clr{3} = {'k-','y-','g-'};
    figure(101);
    for c = 1:length(cohs)
        subplot(spRows,length(cohs),c); hold on;
        for m = 1:length(mods)     % m c d h
            h(m) = plot(hdgs, squeeze(fitInterp.pRight(m,c,D,:)), [clr{c}{m}],'linew',1.5);
            ylim([0 1]);
            text(hdgs(1)+1,1.0-m*0.12,modlabels{m},'color',clr{c}{m}(1),'fontsize',fsz);
        end
        if ~conftask && ~RTtask, xlabel(xLab); end
        if c==1, ylabel('P(Right)'); end
        if length(mods)>1; title(cohlabs{c}); end
   
        set(gca,'xtick',xt);
        set(gca,'ytick',0:0.25:1,'yticklabel',{'0','','0.5','','1'});
        try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

    
        if conftask
            subplot(spRows,length(cohs),c+length(cohs)); hold on;
            for m = 1:length(mods)
                h(m) = plot(hdgs, squeeze(fitInterp.confMean(m,c,D,:)), [clr{c}{m}],'linew',1.5);
                ylim(confYlims); 
            end
            if ~RTtask, xlabel(xLab); end
            if c==1, ylabel(yLab); end
            
            set(gca,'xtick',xt);
            set(gca,'ytick',0:0.25:1,'yticklabel',{'0','0.25','0.5','0.75','1'});
            try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; end

        end
        
        if RTtask
            subplot(spRows,length(cohs),c+length(cohs)*2); hold on;
            for m = 1:length(mods)
                h(m) = plot(hdgs, squeeze(fitInterp.RTmean(m,c,D,:)), [clr{c}{m}],'linew',1.5); 
            end
            xlabel(xLab); 
            if c==1, ylabel('RT (s)'); end
            ylim(RTylims);
            
            set(gca,'xtick',xt);
            try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch;  end

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
            subplot(spRows,length(cohs),c); hold on
            for d = 1:length(deltas)     % m c d h
                h(d) = plot(hdgs, squeeze(fitInterp.pRight(3,c,d,:)), [clr{c}{d}],'linew',1.5); 
                L{d} = sprintf('\x0394=%d',deltas(d));
                text(hdgs(1)+1,1.0-d*0.16,L{d},'color',clr{c}{d}(1),'fontsize',fsz);
            end
            if length(mods)>1; title(cohlabs{c}); end
%             xlabel(xLab); 
            if c==1, ylabel('P(Right)'); end
            ylim([0 1]);
            set(gca,'xtick',xt);
            set(gca,'ytick',0:0.25:1,'yticklabel',{'0','','0.5','','1'});
            try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end


            if conftask>0
                subplot(spRows,length(cohs),c+length(cohs)); hold on
                for d = 1:length(deltas)
                    h(d) = plot(hdgs, squeeze(fitInterp.confMean(3,c,d,:)), [clr{c}{d}],'linew',1.5);
                end
%                 xlabel(xLab); 
                if c==1, ylabel(yLab); end
                ylim(confYlims);
                set(gca,'xtick',xt);
                set(gca,'ytick',0:0.25:1,'yticklabel',{'0','0.25','0.5','0.75','1'});
                try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

            end
            
            
            if RTtask
                subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); hold on;
                for d = 1:length(deltas)
                    h(d) = plot(hdgs, squeeze(fitInterp.RTmean(3,c,d,:)), [clr{c}{d}],'linew',1.5);
                end
                xlabel(xLab); 
                if c==1, ylabel('RT (s)'); end
                ylim(RTylims)
                try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

            end
            
            
        end
        
    end
    
end



