% offline analysis

% human dots3DMP 
% SEP task

clear all; close all

conftask = 1; % 1=colorbars, 2=PDW
normalize = 1;

% specify which (previously saved) mat file to load
 
subject = 'human';
paradigm = 'dots3DMP';
dateRange = 20200213:20210707; % RT
RTtask = 1;

folder = '/Users/stevenjerjian/Desktop/FetschLab/PLDAPS_data/';
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
load([folder file], 'data');

%%
% some new useful vars
for k = 1:length(data.filename)
    data.date(k,1) = str2double(data.filename{k}(9:16));
    data.subjDate{k,:} = [data.subj{k} data.filename{k}(9:16)];
end

%%
% other manual excludes (e.g., RT training)
excludes_filename = {'humanIPQ20200227dots3DMP0904_basic','humanVZC20200229dots3DMP1239_basic'};
excludes_subjDate = {'FRK20200216','FRK20200223','NKT20200215','VZC20200222'};
excludes_subj = {'ASQ', 'HXL', 'XRJ', 'EMF','NEX'};

% excludes_filename = {};
% excludes_subjDate = {};
% excludes_subj = {};

removethese = ismember(data.filename,excludes_filename) | ismember(data.subjDate,excludes_subjDate) | ismember(data.subj,excludes_subj); %#ok<NASGU>
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end
% now this should reflect only good data, per spreadsheet:
blocks = unique(data.filename);

%%
subjs = {'DRH','SJJ','LLV','IPQ','FRK'};
% subjs = {'AAW' 'LLV' 'CXD' 'DRH' 'IPQ' 'SJJ' 'VZC'}; % all 'good' data (pre and post RT)


% remove invalid trials (fixation breaks (which gives nans), excluded subj,
% and obvious testing trials, signaled by very large confidence (saccade
% endpoint) values
removethese = isnan(data.choice) | ~ismember(data.subj,subjs) | data.conf>3 | isnan(data.conf);
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end


% quick look at blocks, for when some need to be excluded
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

% we can be pretty sure blocks with <N trials (say, 30) are to be discarded
removethese = ismember(data.filename,blocks(nTrialsByBlock<30));
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

data.coherence(data.coherence<=0.5) = 0.4;
data.coherence(data.coherence>0.5) = 0.7;
cohs = unique(data.coherence);

% remove the rest
removethese = ~ismember(data.coherence,cohs) & data.modality~=1;
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end
    
% the coh assigned to vestib trials (as a placeholder) depends on which
% cohs were shown in a particular block, so we need to standardize it:
data.coherence(data.modality==1) = cohs(1);

deltas = unique(data.delta); % aka conflict angle
data.heading(abs(data.heading)<0.01) = 0;
hdgs = unique(data.heading);
hdgs(hdgs==0) = [];

% remove the rest
removethese = ~ismember(data.heading,hdgs);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

% final look at blocks 
blocks = unique(data.filename);
nTrialsByBlock = nan(length(blocks),1);
for u = 1:length(blocks)
    nTrialsByBlock(u) = sum(ismember(data.filename,blocks(u)));
end

%% normalize confidence ratings, *within subject*
if normalize

data_orig = data;
usubj = unique(data.subj);
for s = 1:length(usubj)
    data = data_orig;
    removethese = ~strcmp(data.subj,usubj{s});
    for F = 1:length(fnames)
        eval(['data.' fnames{F} '(removethese) = [];']);
    end    
    
    % subtract min and divide by max
    % SJ 07-2021
    % SEPARATELY FOR LEFT AND RIGHT CHOICES 
    L = data.choice==1;
    data.conf(L) = (data.conf(L) - min(data.conf(L))) / max((data.conf(L) - min(data.conf(L))));
    R = data.choice==2;
    data.conf(R) = (data.conf(R) - min(data.conf(R))) / max((data.conf(R) - min(data.conf(R))));

    % OR
    
    % subtract/divide by *means*
%     minOfMeans = nanmin(nanmin(nanmin(nanmin(confMean))));
%     data.conf = data.conf - minOfMeans;
%     dots3DMP_parseData
%     maxOfMeans = nanmax(nanmax(nanmax(nanmax(confMean))));
%     data.conf = data.conf / maxOfMeans;

    % OR
    
    % simply cap at 1/0
    % data.conf(data.conf>1) = 1;
    % data.conf(data.conf<0) = 0;

    % append each subj to a new data struct
    if s==1
        data_new = data;
    else
        for F = 1:length(fnames)
            eval(['data_new.' fnames{F} '(end+1:end+length(data.date)) = data.' fnames{F} ';']);
        end
    end
end
data = data_new;
clear data_new data_orig
end


%%

save([file(1:end-4) '_clean.mat'],'data')