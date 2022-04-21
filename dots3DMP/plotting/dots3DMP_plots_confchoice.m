function dots3DMP_plots_confchoice(parsedData,mods,cohs,deltas,conftask)

if conftask == 0
    disp('Cannot plot conf vs choice for non-conf task, doh!')
    return
end

fsz = 15; % fontsize

% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
D = find(deltas==0);

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
% clr{1} = {'ko','ro','bo'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(111);
% set(gcf,'Color',[1 1 1],'Position',[300 1000 450+300*(length(cohs)-2) 200+150*(conftask>0)+150*RTtask],'PaperPositionMode','auto'); clf;
set(gcf,'Color',[1 1 1],'Position',[200 80 600 250],'PaperPositionMode','auto'); clf;
subplot(121);
for m = 1:length(mods)
    for c = 1:length(cohs)
    h(m) = plot(squeeze(parsedData.pRight(m,c,D,:)), squeeze(parsedData.confMean(m,c,D,:)), [clr{c}{m} '-'],'linewidth',1.5); hold on;
    end
end
ylim([0.3 1]);
xlabel('P(right)');
if conftask==1, ylabel('SaccEP');
elseif conftask==2, ylabel('P(High Bet)');
end
try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

% now do deltas

if length(deltas)>1
    
    subplot(122);
    
    clr{1} = {'bs','cs','gs'};
    clr{2} = {'b^','c^','g^'};
    clr{3} = {'bo','co','go'};
    linstl = '-:';
    
    clear L;
    for d = 1:length(deltas)
        for c = 1:length(cohs)
            
            hd(c) = plot(squeeze(parsedData.pRight(3,c,d,:)), squeeze(parsedData.confMean(m,c,d,:)), [clr{c}{d} linstl(c)],'linewidth',1.5); hold on;
            
        end
    end
    
    ylim([0.2 1]);
    xlabel('P(right)');
    if conftask==1, ylabel('SaccEP');
    elseif conftask==2, ylabel('P(High Bet)');
    end
    try changeAxesFontSize(gca,fsz,fsz); tidyaxes(gca,fsz); catch; disp('plot clean up skipped'); end

end