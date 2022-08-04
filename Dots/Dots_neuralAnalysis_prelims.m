
% for tuning functions
tuning_vonMises = @(b,dir) b(1) * exp(b(2)*cosd(dir-b(3))) / (2*pi*besseli(0,b(2))) + b(4);
tuning_vonMises_err = @(b,dir,FR) nansum((tuning_vonMises(b,dir)-FR).^2); % sum squared error

% check max trial len, nunits, ntrials (for preallocation)
for n = 1:length(dataCell)
    if max(dataCell{n}.Exp.openEvents.choice)==2 % convert choices to 0:1
        dataCell{n}.Exp.openEvents.choice = dataCell{n}.Exp.openEvents.choice-1;
    end
    if any(isnan(dataCell{n}.Exp.openEvents.choice))
        keyboard
    end
    if any(isnan(dataCell{n}.Exp.openEvents.pdw))
        keyboard
    end
    
    mStart = dataCell{n}.Mapping.openEvents.motionStart';
    mEnd = dataCell{n}.Mapping.openEvents.motionEnd';
    dur = round((mEnd-mStart)*1000);
% % %     maxlen_map(n) = round(max(dur));
    pct99_map(n) = ceil(prctile(dur,99)); % use this instead
    if pct99_map(n)>600; pct99_map(n) = ceil(prctile(dur,98)); end  % kluge one weird session
    
    mStart = dataCell{n}.Exp.openEvents.motionStart';
    mEnd = dataCell{n}.Exp.openEvents.motionEnd';
    dur = round((mEnd-mStart)*1000);
    maxlen_exp(n) = round(max(dur));
    pct99_exp(n) = ceil(prctile(dur,99)); % use this instead

    nU(n) = length(dataCell{n}.Mapping.spikeTimes);
    nTrials_map(n) = length(dataCell{n}.Mapping.openEvents.motionStart);
    nTrials_exp(n) = length(dataCell{n}.Exp.openEvents.motionStart);
end
% % this tells me that 1600 + sum(extRaster) is enough for exp,
% % and 400 + sum(extRaster) is enough for mapping. set these manually and
% % exclude trials that exceed them.
% % the max nch and ntr will be set below.
% rasterLen_exp = 1600;
% rasterLen_map = 400;

rasterLen_exp = max(pct99_exp)+1; % +1 for a cushion, some weird edge bugs
rasterLen_map = max(pct99_map)+1;

nUnit = sum(nU) - size(excludes,1);
sessID = nan(nUnit,1); % keep track of which session a given unit was in

% preallocate

% raster: [unit,trials,time]
raster_map =    nan(nUnit,max(nTrials_map),rasterLen_map+sum(extRaster));
raster_dotsOn = nan(nUnit,max(nTrials_exp),rasterLen_exp+sum(extRaster));
raster_RT =     nan(nUnit,max(nTrials_exp),rasterLen_exp+sum(extRaster));

% aggregated, normalized psth [unit,conditions(x4),time]:
psthNorm_dotsOn_one =   nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_two =   nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_three = nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_dotsOn_four =  nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_one =       nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_two =       nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_three =     nan(nUnit,rasterLen_exp+sum(extRaster));
psthNorm_RT_four =      nan(nUnit,rasterLen_exp+sum(extRaster));

% spike rates/counts: [unit,trials]
% % % spRate_map =  nan(nUnit,max(nTrials_map));
spRate_exp =  nan(nUnit,max(nTrials_exp));
spCount_exp = nan(nUnit,max(nTrials_exp));

% misc
prefDir_fit = nan(nUnit,1);
CP_high =     nan(nUnit,1);
CP_low =      nan(nUnit,1);
CP_all =      nan(nUnit,1);
ConfP_pref =  nan(nUnit,1);
ConfP_null =  nan(nUnit,1);
ConfP_all =   nan(nUnit,1);

