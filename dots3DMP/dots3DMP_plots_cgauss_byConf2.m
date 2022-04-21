function dots3DMP_plots_cgauss_byConf2(gfit,parsedData,mods,cohs,deltas,hdgs,conftask,RTtask)

% loop over 'conditions' indexed as [modality coherence], with the
% familiar manual kluge necessary to avoid invalid combinations:
ucoh = unique(data.coherence);
ucond = [1 ucoh(1); 2 ucoh(1); 2 ucoh(2); 3 ucoh(1); 3 ucoh(2)];
titles = {'Ves';'Vis-Low';'Vis-High';'Comb-Low';'Comb-High';'All'};
subplotInd = [2 3 4 5 6 1];

D = length(deltas)+1; % (the extra column we made for pooling across deltas)

% CHOICES i.e. pRight
figure(201);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c=1:length(ucond)+1
    subplot(3,2,subplotInd(c));

    for cc=1:2
        beta = [gfit.choice.mu(m,c,D,cc) gfit.choice.sigma(m,c,D,cc)];
        h(cc) = plot(parsedData.xVals, gfit.choice.func(beta,parsedData.xVals), clr{cc}{m}(1),'linewidth',1.5); hold on;
        errorbar(hdgs, squeeze(parsedData.pRight(m,c,D,:,cc)), squeeze(parsedData.pRightSE(m,c,D,:,cc)), clr{cc}{m},'linewidth',1.5);
    end
           
        ylim([0 1]);
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
        else
            title(modlabels{m})
        end
    legend(h,'High Bet','Low Bet','Location','northwest');
    if m==length(mods), xlabel(xLab); end
    if m==2, ylabel('P(right)'); end
    try changeAxesFontSize(gca,15,15); catch; end
    end
end

% RT

if RTtask
    
figure(201+D+1);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
%         subplot(length(mods),length(cohs),c+(m-1)*length(cohs))
        subplot(1,length(mods),m)
        if m==1 & c<1,continue,end
        
        for cc=1:2
            beta = [gfit.RT.ampl(m,c,D,cc) gfit.RT.mu(m,c,D,cc) gfit.RT.sigma(m,c,D,cc) gfit.RT.bsln(m,c,D,cc)];
            h(cc) = plot(parsedData.xVals, gfit.RT.func(beta,parsedData.xVals), clr{cc}{m}(1),'linewidth',1.5); hold on;
            errorbar(hdgs, squeeze(parsedData.RTmean(m,c,D,:,cc)), squeeze(parsedData.RTse(m,c,D,:,cc)), clr{cc}{m},'linewidth',1.5);
        end
           
        ylim([0.3 1]);
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
        else
            title(modlabels{m})
        end
    legend(h,'High Bet','Low Bet','Location','northwest');
    if m==length(mods), xlabel(xLab); end
    if m==2, ylabel('RT (s)'); end
    try changeAxesFontSize(gca,15,15); catch; end
    end
end

end

%{
% CONFIDENCE (SANITY CHECK)
figure(201+D+2);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c=1:length(cohs)
    for m=1:length(mods)
        subplot(length(mods),length(cohs),c+(m-1)*length(cohs)); hold on
        if m==1 && c~=1, delete(gca); continue, end
        
        for cc=1:2
            errorbar(hdgs, squeeze(parsedData.confMean(m,c,D,:,cc)), squeeze(parsedData.confSE(m,c,D,:,cc)), clr{c}{m},'linewidth',1.5);
        end
           
        ylim([0 1]);
        if m~=1
            title([modlabels{m} ', coh = ' num2str(cohs(c))]); 
        else
            title(modlabels{m})
        end
    if m==length(mods), xlabel('heading angle (deg)'); end
    if m==2
        if conftask==1, ylabel('P(High Bet'); 
        else,           ylabel('SEP')
        end
    end
    try changeAxesFontSize(gca,15,15); catch; end
    end
end
%}