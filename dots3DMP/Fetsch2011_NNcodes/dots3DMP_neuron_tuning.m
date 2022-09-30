function [meanFRs,semFRs] = dots3DMP_neuron_tuning(data,mods,cohs,deltas,hdgs)

numunits = numel(unique(data.unitnum));
n = nan(length(mods),length(cohs),length(deltas)+1,length(hdgs),numunits);

meanFRs = n;
semFRs  = n;

for u = 1:numunits
for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)+1 % add extra column for all trials irrespective of delta
    for h = 1:length(hdgs)
        if d==length(deltas)+1
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h); % all trials irrespective of delta
        else
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        end
        
        % get data just for one unit
        J = J & (data.unitnum == u);
        
        n(m,c,d,h,u) = nansum(J);
        %pRight(m,c,d,h) = nansum(J & data.choice==1) / n(m,c,d,h); % 1 is rightward!!
        fr = data.spikeRates(J);
        if isempty(fr), continue, end

        meanFRs(m,c,d,h,u) = mean(fr);
        semFRs(m,c,d,h,u) = std(fr) ./ sqrt(n(m,c,d,h,u));
    end
end
end
end
end

for c = 1:length(cohs)
    n(1,c,:,:,:) = n(1,1,:,:,:);
    meanFRs(1,c,:,:,:) = meanFRs(1,1,:,:,:);
    semFRs(1,c,:,:,:) = semFRs(1,1,:,:,:);
end
    
