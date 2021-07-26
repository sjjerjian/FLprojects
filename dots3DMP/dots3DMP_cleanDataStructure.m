function [data, nTrialsByDate, nTrialsByBlock] = dots3DMP_cleanDataStructure(data,excludes_filename,excludes_date,minBlockTrialCount)

% clean up the data structure to remove 

% some new useful vars
for k = 1:length(data.filename)
    data.date(k,1) = str2double(data.filename{k}(6:13));
    data.subjDate{k,:} = [data.subj{k} data.filename{k}(6:13)];
end


%% Some manual excludes e.g. of bad sessions, non-RT/PDW

% excludes_filename = {};
% excludes_date = [20210608];

% remove fixation breaks (indicated by nans), or one-target trials
% also remove excludes based on filenames and dates
removethese = isnan(data.choice) | isnan(data.RT) | isinf(data.RT) | isnan(data.PDW) | data.oneConfTargTrial | data.oneTargTrial;
removethese = removethese | ismember(data.filename,excludes_filename) | ismember(data.date,excludes_date);
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% should do the trial number based exclusion here, earlier on we are
% counting fixation breaks

% quick look at blocks, for when some need to be excluded
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

% we can be pretty sure blocks with <N good trials are to be discarded
% N = 50;
removethese = ismember(data.filename,blocks(nTrialsByBlock<minBlockTrialCount));
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end
% quick look at blocks again
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end


%% cull data

mods = unique(data.modality); 

data.delta(data.delta==-2) = -3;
data.delta(data.delta==2) = 3;
deltas = unique(data.delta); % aka conflict angle

% simplify cohs (collapse similar ones)
data.coherence(data.coherence<=0.5) = 0.4;
data.coherence(data.coherence>0.5) = 0.8;
cohs = unique(data.coherence);

% the coh assigned to vestib trials (as a placeholder) depends on which
% cohs were shown in a particular block, so we need to standardize it:
data.coherence(data.modality==1) = cohs(1);

data.heading(abs(data.heading)<0.01) = 0;
% hdgs = unique(data.heading);
hdgs = [-12 -6 -3 -1.5 0 1.5 3 6 12];

% remove the rest
removethese = ~ismember(data.heading,hdgs) | ~ismember(data.coherence,cohs) | ~ismember(data.delta,deltas);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% remove unnecessary fields and order fields
data = rmfield(data,{'reward','subj','oneTargTrial','oneConfTargTrial','TargMissed','subjDate'});

sorted_fnames = {'filename','date','heading','modality','coherence','delta','choice','RT','PDW','correct'};
data = orderfields(data,sorted_fnames);


%% check block and day trial counts
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

dates = unique(data.date);
nTrialsByDate = nan(length(dates),1);
for u = 1:length(dates)
    nTrialsByDate(u) = sum(ismember(data.date,dates(u)));
end




