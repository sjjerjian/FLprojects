% confidence (PDW) for correct and error trials as a function of heading angle
% ignore other stim conditions for now

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
xlabel('heading (deg)');
changeAxesFontSize(gca,14,14);

if conftask==1
    ylabel('Sacc EndPoint');
    ylim([-0.3 1.3])
elseif conftask==2
    ylabel('P(HighBet)');
    ylim([0.2 0.9])
end
box off;

clear confCorr confErr nCorr nErr
ushdgs = hdgs(hdgs>0);
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
    
    RTCorr(h) = mean(data.RT(I));
    RTCorrSE(h) = std(data.RT(I))/sqrt(sum(I));
    
    RTErr(h) = mean(data.RT(J));
    RTErrSE(h) = std(data.RT(J))/sqrt(sum(J));
end

if conftask==2
confCorrSE = sqrt( (confCorr.*(1-confCorr)) ./ nCorr);
confErrSE = sqrt( (confErr.*(1-confErr)) ./ nErr);
end

subplot(132); plot(ushdgs,confCorr,'b-o',ushdgs,confErr,'r-o','linew',1.5);
xlabel('|heading| (deg)');
changeAxesFontSize(gca,14,14); box off;
text(6,0.7,'Correct','color','b','fontsize',16);
text(6,0.64,'Error','color','r','fontsize',16);
for h=1:length(ushdgs)
    text(ushdgs(h),confCorr(h)+0.02,num2str(nCorr(h)),'fontsize',12,'horizo','center');
    text(ushdgs(h),confErr(h)+0.02,num2str(nErr(h)),'fontsize',12,'horizo','center');
end

subplot(133); 
plot(RTCorr,confCorr,'b-o',RTErr,confErr,'r-o','linew',1.5);
xlabel('RT');
changeAxesFontSize(gca,14,14); box off;
for h=1:length(ushdgs)
    text(RTCorr(h)-0.003,confCorr(h)-0.05,sprintf('%.1f',ushdgs(h)),'fontsize',12,'horizo','center');
    text(RTErr(h)-0.003,confErr(h)-0.05,sprintf('%.1f',ushdgs(h)),'fontsize',12,'horizo','center');
end
if conftask==1
    xlim([0.7 1.5])
elseif conftask==2
    xlim([0.57 0.74])
end

%% Confidence for different deltas as a function of heading angle

for h = 1:length(hdgs)
    for d=1:length(deltas)
        I = data.heading==hdgs(h) & data.modality==3 & data.delta==deltas(d);
        nTr(h,d) = sum(I);
    
        if conftask==1
            confdelta(h,d) = mean(data.conf(I));
            confdeltaSE(h,d) = std(data.conf(I))/sqrt(sum(I));
        elseif conftask==2
            confdelta(h,d) = mean(data.PDW(I));
        end
    end
end

cols = {'b','c','g'};
figure('color','white','position',[300 300 700 250]);
subplot(131); hold on;
for d=1:length(deltas)
    plot(hdgs,confdelta(:,d),'linew',1.5,'marker','o','color',cols{d});
end
xlabel('heading (deg)');
changeAxesFontSize(gca,14,14);

if conftask==1
    ylabel('Sacc EndPoint');
    ylim([-0.3 1.3])
elseif conftask==2
    ylabel('P(HighBet)');
    ylim([0.3 1.0])
end
box off;

clear confdelta* RTdelta* nTr
ushdgs = hdgs(hdgs>0);
for h = 1:length(ushdgs)   
    for d=1:length(deltas)
    I = abs(data.heading)==ushdgs(h) & data.modality==3 & data.delta==deltas(d);
    nTr(h,d) = sum(I);
    
    if conftask==1
        confdelta(h,d) = mean(data.conf(I));
        confdeltaSE(h,d) = std(data.conf(I))/sqrt(sum(I));
    elseif conftask==2
        confdelta(h,d) = mean(data.PDW(I));
    end

    RTdelta(h,d) = mean(data.RT(I));
    RTdeltaSE(h,d) = std(data.RT(I))/sqrt(sum(I));
    end
end

if conftask==2
confdeltaSE = sqrt( (confdelta.*(1-confdelta)) ./ nTr);
end

subplot(132); hold on;
for d=1:length(deltas)
    plot(ushdgs,confdelta(:,d),'linew',1.5,'marker','o','color',cols{d});
%     for h=1:length(ushdgs)
%         text(ushdgs(h),confdelta(d)+0.02,num2str(nTr(h,d)),'fontsize',12,'horizo','center');
%     end
end
xlabel('|heading| (deg)');
changeAxesFontSize(gca,14,14); box off;
ylim([0.3 1])
% text(6,0.7,'Correct','color','b','fontsize',16);
% text(6,0.64,'Error','color','r','fontsize',16);


subplot(133); hold on;
for d=1:length(deltas)
    plot(RTdelta(:,d),confdelta(:,d),'linew',1.5,'marker','o','color',cols{d});
%     for h=1:length(ushdgs)
%         text(RTdelta(h,d)-0.003,confdelta(h,d)-0.05,sprintf('%.1f',ushdgs(h)),'fontsize',12,'horizo','center');
%     end
%     
end
xlabel('RT');
if conftask==1
    xlim([0.7 1.5])
elseif conftask==2
    xlim([0.57 0.74])
end
ylim([0.3 1])
changeAxesFontSize(gca,14,14); box off;
