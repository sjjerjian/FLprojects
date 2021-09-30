function dots3DMP_CorrectVsErrorCurves(data,hdgs,mods,cohs,conftask,RTtask)
% look at PDW and RT separately for correct and incorrect trials

removethese = ~ismember(data.modality,mods) | ~ismember(data.coherence,cohs);
% removethese = removethese | data.heading==0 | abs(data.delta)>0;

fnames = fieldnames(data);
for f=1:length(fnames)
    data.(fnames{f})(removethese) = [];
end

for h = 1:length(hdgs)
    I = data.heading==hdgs(h) & data.correct;
    nCorr(h) = sum(I);
    
    if conftask==1
        confCorr(h) = mean(data.conf(I));
        coffCorrSE(h) = std(data.conf(I))/sqrt(sum(I));
    elseif conftask==2
        confCorr(h) = mean(data.PDW(I));
    end
    
    J = data.heading==hdgs(h) & ~data.correct;
    nErr(h) = sum(J);
    if conftask==1
        confErr(h) = mean(data.conf(J));
        coffErrSE(h) = std(data.conf(J))/sqrt(sum(J));
    elseif conftask==2       
        confErr(h) = mean(data.PDW(J));
    end
end

figure('color','white','position',[300 300 700 250]);
subplot(131); plot(hdgs,confCorr,'b-o',hdgs,confErr,'r-o','linew',1.5);
xlabel('heading (deg)'); xtickangle(45);
set(gca,'xtick',hdgs,'tickdir','out');
try changeAxesFontSize(gca,12,12); offsetAxes; catch; end

ydata = [confCorr confErr];
if conftask==1
    ylabel('SEP');
    ylim([-0.3 1.3])
elseif conftask==2
    ylabel('P(HighBet)');
    ylim([0.9 1.1].*([min(ydata(:)) max(ydata(:))]))
end
box off;

clear confCorr confErr nCorr nErr
ushdgs = hdgs(hdgs>=0);
for h = 1:length(ushdgs)   
    I = abs(data.heading)==ushdgs(h) & data.correct;
    nCorr(h) = sum(I);
    
    if conftask==1
        confCorr(h) = mean(data.conf(I));
        coffCorrSE(h) = std(data.conf(I))/sqrt(sum(I));
    elseif conftask==2
        confCorr(h) = mean(data.PDW(I));
    end
     
    J = abs(data.heading)==ushdgs(h) & ~data.correct;
    nErr(h) = sum(J);
    if conftask==1
        confErr(h) = mean(data.conf(J));
        coffErrSE(h) = std(data.conf(J))/sqrt(sum(J));
    elseif conftask==2
        confErr(h) = mean(data.PDW(J));
    end
    
    if RTtask
        RTCorr(h) = mean(data.RT(I));
        RTCorrSE(h) = std(data.RT(I))/sqrt(sum(I));
        
        RTErr(h) = mean(data.RT(J));
        RTErrSE(h) = std(data.RT(J))/sqrt(sum(J));
    end
end

if conftask==2
    confCorrSE = sqrt( (confCorr.*(1-confCorr)) ./ nCorr);
    confErrSE = sqrt( (confErr.*(1-confErr)) ./ nErr);
end


subplot(1,2+double(conftask>0),2); plot(ushdgs,confCorr,'b-o',ushdgs,confErr,'r-o','linew',1.5);
set(gca,'xtick',ushdgs,'tickdir','out');
xlabel('|heading| (deg)'); xtickangle(45);
try changeAxesFontSize(gca,12,12); offsetAxes; catch; end
box off;
text(6,0.7,'Correct','color','b','fontsize',14);
text(6,0.64,'Incorrect','color','r','fontsize',14);
for h=1:length(ushdgs)
    text(ushdgs(h),confCorr(h)+0.02,num2str(nCorr(h)),'fontsize',12,'horizo','center');
    text(ushdgs(h),confErr(h)+0.02,num2str(nErr(h)),'fontsize',12,'horizo','center');
end
ydata = [confCorr confErr];
ylim([0.9 1.1].*([min(ydata(:)) max(ydata(:))]))


if RTtask
    subplot(1,3,3);
    plot(RTCorr,confCorr,'b-o',RTErr,confErr,'r-o','linew',1.5);
    xlabel('RT');
    xlim([.5 .8]); ylim([0.9 1.1].*([min(ydata(:)) max(ydata(:))]))
    set(gca,'xtick',0.5:0.05:0.8,'xticklabel',{'.5','','.6','','.7','','.8'},'tickdir','out');
    try changeAxesFontSize(gca,12,12); offsetAxes; catch; end
    box off;
    for h=1:length(ushdgs)
        text(RTCorr(h)+0.03,confCorr(h),sprintf('%g',ushdgs(h)),'fontsize',12,'horizo','center');
        text(RTErr(h)+0.03,confErr(h),sprintf('%g',ushdgs(h)),'fontsize',12,'horizo','center');
    end
end


% if conftask==1
%     xlim([0.7 1.5])
% elseif conftask==2
%     xlim([0.57 0.74])
% end