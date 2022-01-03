function dots3DMP_plots_cgauss_byCoh(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)
% SJ 07-2021 converted to function for cleaner workspace

% first, for all trials irrespective of delta
% D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);
D = gfit.D;

if all(mods==1), cohs=1; end
fsz = 15; % fontsize

spRows = 1 + double(conftask>0) + double(RTtask);

modlabels = {'Ves','Vis','Comb'};
xLab = sprintf('heading angle (%s)',char(176));

if conftask==1 
    yLab = 'SaccEP'; confYlims = [0.2 0.9]; RTyt = 0:0.25:2; 
    if all(mods==1), RTylims = [0.9 1.5]; 
    else,            RTylims = [0.9 1.7]; 
    end
    xt = -10:5:10;
     if length(mods)>1, cohlabs = {'Low Coh','High Coh'}; end
elseif conftask==2
    yLab = 'P(High Bet)'; confYlims = [0.4 1.0]; RTyt = 0:0.1:2; 
    if all(mods==1), RTylims = [0.5 0.72]; 
    else,            RTylims = [0.5 0.9]; 
    end
    xt = -12:6:12;
    if length(mods)>1, cohlabs = {sprintf('coh = %.1f',cohs(1)),sprintf('coh = %.1f',cohs(2))}; end
end 


         %ves %vis %comb
% clr{1} = {'ko','mo','co'};
clr{1} = {'ko','ro','bo'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
clr{4} = clr{1};
figure(101+D);
% set(gcf,'Color',[1 1 1],'Position',[300 1000 230+300*(length(cohs)-1) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
set(gcf,'Color',[1 1 1],'Position',[200 80 700 900],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    % choice
    subplot(spRows,length(cohs),c); box off; hold on;
    for m = 1:length(mods)     % m c d h
        beta = [gfit.choice.mu(m,c,D) gfit.choice.sigma(m,c,D)];
        h(m) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.pRightSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        text(hdgs(1)+1,1.0-m*0.12,modlabels{m},'color',clr{c}{m}(1),'fontsize',fsz);
%         text(hdgs(1)+0.5,1.0-m*0.07,sprintf('%s: mu = %.2f, s = %.2f',modlabels{m},beta(1),beta(2)),'color',clr{c}{m}(1))
    end
    if length(mods)>1; title(cohlabs{c}); end
    ylim([0 1]);
    set(gca,'xtick',xt);
    set(gca,'ytick',0:0.25:1,'yticklabel',{'0','.25','.5','.75','1'});
    if ~conftask && ~RTtask, xlabel(xLab); end
    if c==1, ylabel('P(right)'); end
    try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end
    
        % conf
    if conftask
        subplot(spRows,length(cohs),c+length(cohs)); box off; hold on;
        for m = 1:length(mods)
            beta = [gfit.conf.ampl(m,c,D) gfit.conf.mu(m,c,D) gfit.conf.sigma(m,c,D) gfit.conf.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:)), squeeze(parsedData.confSE(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        end
%         if length(mods)>1; title(cohlabs{1}); end
        set(gca,'xtick',xt);
        set(gca,'ytick',0:0.25:1,'yticklabel',{'0','.25','.5','.75','1'});

%         ylim(confYlims);
%         xlabel(xLab); 
        if c==1, ylabel(yLab); end

        try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end
    end
    
    % RT
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); box off; hold on;
        for m = 1:length(mods)        
            beta = [gfit.RT.ampl(m,c,D) gfit.RT.mu(m,c,D) gfit.RT.sigma(m,c,D) gfit.RT.bsln(m,c,D)];
            h(m) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{c}{m}(1) '-'],'linewidth',1.5); hold on;       
            errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:)), squeeze(parsedData.RTse(m,c,D,:)), clr{c}{m},'linewidth',1.5);
        end
%         if length(mods)>1; title(cohlabs{1}); end
%         set(gca,'xtick',xt,'ytick',0:.1:2);
%         ylim(RTylims)
        xlabel(xLab); 
        if c==1, ylabel('RT (s)'); end
        try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end
    end
end


%% now separate by delta

if length(deltas)>1 && any(mods==3)
    
clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};
clr{4} = {'bd','cd','gd'};


% clr{1} = {'bs','gs'};
% clr{2} = {'b^','g^'};
% clr{3} = {'bo','go'};
% clr{4} = {'bd','gd'};


clear L;
figure(208);
% set(gcf,'Color',[1 1 1],'Position',[50 20 230+300*(length(cohs)-1) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
set(gcf,'Color',[1 1 1],'Position',[900 80 700 900],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    % choice
    subplot(spRows,length(cohs),c); box off; hold on;
    for d = 1:length(deltas)     % m c d h
        beta = [gfit.choice.mu(3,c,d) gfit.choice.sigma(3,c,d)];
        h(d) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.pRightSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
        L{d} = sprintf('\\Delta=%d',deltas(d));
        text(hdgs(1)+1,1.0-d*0.16,L{d},'color',clr{c}{d}(1),'fontsize',fsz);
    end
    if length(mods)>1; title(cohlabs{c}); end
    set(gca,'xtick',xt);
    set(gca,'ytick',0:0.25:1,'yticklabel',{'0','.25','.5','.75','1'});
    ylim([0 1]);
%     lh=legend(h,L,'location','southeast'); set(lh,'box','off');
%     xlabel(xLab); 
    if c==1, ylabel('P(Right)'); end
    try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

    % conf
    if conftask
        subplot(spRows,length(cohs),c+length(cohs)); box off; hold on;
        for d = 1:length(deltas)
            beta = [gfit.conf.ampl(3,c,d) gfit.conf.mu(3,c,d) gfit.conf.sigma(3,c,d) gfit.conf.bsln(3,c,d)];
            h(d) = plot(parsedData.xVals, gfit.conf.func(beta,parsedData.xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;       
            errorbar(hdgs, squeeze(parsedData.confMean(3,c,d,:)), squeeze(parsedData.confSE(3,c,d,:)), clr{c}{d},'linewidth',1.5);
            L{d} = sprintf('?=%d',deltas(d));
        end
%         if length(mods)>1; title(cohlabs{c}); end
        set(gca,'xtick',xt);
%         ylim(confYlims);
    %     legend(h,L,'location','northwest');
%         xlabel(xLab); 
        if c==1, ylabel(yLab); end
        set(gca,'ytick',0:0.25:1,'yticklabel',{'0','.25','.5','.75','1'});
        try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end
    end

    % RT
    if RTtask
        subplot(spRows,length(cohs),c+length(cohs)*(2-(conftask==0))); box off; hold on;
        for d = 1:length(deltas)
            beta = [gfit.RT.ampl(3,c,d) gfit.RT.mu(3,c,d) gfit.RT.sigma(3,c,d) gfit.RT.bsln(3,c,d)];
            h(d) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), [clr{c}{d}(1) '-'],'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.RTmean(3,c,d,:)), squeeze(parsedData.RTse(3,c,d,:)), clr{c}{d},'linewidth',1.5);
            L{d} = sprintf('?=%d',deltas(d));
        end
%         if length(mods)>1; title(cohlabs{c}); end
%         set(gca,'xtick',xt,'ytick',0:.1:2);
        xlabel(xLab); 
        if c==1, ylabel('RT (s)'); end
%         ylim(RTylims);
        try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz);catch; end
    end    

end

% suptitle('Combined condition: Cue Conflict')

end