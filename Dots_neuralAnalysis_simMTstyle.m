% 'sim' style decode of real MT data



%% convert to format for simMT

T = 1; U = 0;
raster = nan(40691,size(raster_dotsOn,4)-extRaster(1));
for n = 1:length(dataCell)
    n
    for c = 1:nCh(n)
        U = U+1;
        for t = 1:length(dataCell{n}.Exp.openEvents.coherence)
            if ~isnan(spRate_exp(n,c,t)) && ~isnan(spCount_exp(n,c,t))
%                 sprate(T,1) = spRate_exp(n,c,t); 
%                 spcount(T,1) = spCount_exp(n,c,t);
                raster(T,:) = raster_dotsOn(n,c,t,extRaster(1)+1:end);

                unit(T,1) = U; 
                cohD(T,1) = dataCell{n}.Exp.openEvents.coherence(t);
                dirD(T,1) = dataCell{n}.Exp.openEvents.direction(t); 
                choiceD(T,1) = dataCell{n}.Exp.openEvents.choice(t); 
                pdwD(T,1) = dataCell{n}.Exp.openEvents.pdw(t);            
                mStart = dataCell{n}.Exp.openEvents.motionStart(t);            
                mEnd = dataCell{n}.Exp.openEvents.motionEnd(t);
                dur(T,1) = round((mEnd-mStart)*1000);        
                T = T+1;
            end
        end   
    end
end

% this kind of pseudopop approach relies on sampling from the neural
% responses for a set of N simulated trials.

% what I'm missing is a way to normalize spike trains, but let's try
% without it for now [can normalize PSTHs, that should work]


%% 'build'

nTrials = 1000;

% signed motion coherence; negative is leftward
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,nTrials,'true')';
uscoh = abs(coh);
dir = nan(size(coh));
dir(coh<0) = 180;
dir(coh>0) = 0;
dir(coh==0) = randsample([0 180],sum(coh==0),'true');

% R = [trial, neuron, time]
R = nan(nTrials,max(unit),size(raster,2));
for t = 1:nTrials
    for U = 1:max(unit)
        matchingTrials = find(unit==U & cohD==abs(coh(t)) & dirD==dir(t));
        draw = randsample(matchingTrials,1); % a random sample for each cell, so running this many times will bootstrap any result
        R(t,U,:) = raster(draw,:);
    end
end


%% 'decode'

latencyOrig = latency;
latency = latency(1);
prefDirs = prefDir_fit(~isnan(prefDir_fit));


simMT_decode_DDM


