

% need to bring in nexonar Datapixx timestamps from PLDAPS file

ntrs = length(nex.nexdata);
hdgs = unique(nex.conditions.heading);

mcols = {'Greys','Reds','Blues'};

m = 1;
cmap = flipud(cbrewer('seq',mcols{m},length(hdgs)*3));
cmap = cmap((1:length(hdgs))*2,:);

trcount=0;
figure; hold on;
for tr = 1 : ntrs-1
    
    if nex.behavior.goodtrial(tr)==0 || nex.conditions.modality(tr)~=m, continue, end
    
    h = hdgs == nex.conditions.heading(tr) & nex.conditions.delta(tr)==0;
    if sum(h)
        trcount=trcount+1;
        plot(nex.nexdata{tr}(:,1),nex.nexdata{tr}(:,3)-mean(nex.nexdata{tr}(1:50,3)),'color',cmap(h,:),'linew',2)
%         plot(gradient(nex.nexdata{tr}(:,3)-mean(nex.nexdata{tr}(1:50,3))),'color',cmap(h,:),'linew',2)
    end
end
% ylim([-1 1]*40)