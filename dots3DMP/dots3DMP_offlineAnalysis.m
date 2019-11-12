% offline analysis wrapper for PLDAPS data, dots3DMP paradigm
% CF started it 10-3-18

clear all; close all

conftask=1; % 1=colorbars, 2=PDW

% these will specify which (previously saved) mat file to load
subject = 'human';
dateRange = 20190612:20191231; % everything
folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];

load([folder file], 'data');

% this could be useful
for k = 1:length(data.filename)
    data.date(k,1) = str2double(data.filename{k}(9:16));
end

% struct data has fields:
% filename
% subj: subject code
% choice: 1=left, 2=right, nan = fixation break or otherwise invalid
% heading: angle in deg relative to straight ahead, positive = rightward
% coherence: dot motion coherence aka visual reliability
% modality: stimulus modality: 1=ves, 2=vis, 3=comb
% delta: conflict angle in deg, positive means visual right, vestib left
% correct: was the choice correct, 1 or 0
% conf: confidence rating via saccadic end-point to color-bar target
%       (or in other datasets, PDW)


%% choose subjects to include based on 3-letter code

subjs = unique(data.subj); % all
% subjs = {'AAW'}; % the best single subject
% subjs = {'LLV'};
% subjs = {'AAW' 'LLV'};
% subjs = {'AAW' 'LLV' 'IWT'}; % these other two have the most trials after AAW
% subjs = {'AAW' 'LLV' 'IWT' 'CXD'}; % including CXD would add a few extra trials
% subjs = {'AAW' 'LLV' 'IWT' 'EMF'}; % including CXD would add a few extra trials
% subjs = {'AAW' 'LLV' 'IWT' 'CXD' 'EMF'}; % including CXD would add a few extra trials


% remove invalid trials (fixation breaks (which gives nans), excluded subj,
% and obvious testing trials, signaled by very large confidence (saccade
% endpoint) values
removethese = isnan(data.choice) | ~ismember(data.subj,subjs) | data.conf>2 | isnan(data.conf);
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end


%% parse data

mods = unique(data.modality); 

% [C,IA,IC] = unique(data.coherence)
% N = hist(data.coherence,unique(data.coherence))

% % % cohs = unique(data.coherence); 
    % currently too many cohs in the dataset, so...

% N = hist(data.coherence,unique(data.coherence))
% unique(data.coherence)

% ...lump coherences together (work in progress)

% avery
data.coherence(data.coherence<0.2) = 0.1;
data.coherence(data.coherence>=0.2 & data.coherence<0.6) = 0.5;
data.coherence(data.coherence>=0.6) = 0.9;
cohs = [0.1 0.5 0.9];

% % everyone else
% data.coherence(data.coherence<=0.3) = 0.1;
% data.coherence(data.coherence>0.3 & data.coherence<=0.7) = 0.5;
% cohs = [0.1 0.5];



% remove the rest
removethese = ~ismember(data.coherence,cohs) & data.modality~=1;
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end
    
% the coh assigned to vestib trials (as a placeholder) depends on which
% cohs were shown in a particular block, so we need to standardize it:
data.coherence(data.modality==1) = cohs(1);

    
deltas = unique(data.delta); % aka conflict angle
% % % hdgs = unique(data.heading);
      % same here.  the 1.5-12 range was only used rarely, and in fact is a
      % good signature of warmup or testing-mode trials to be excluded
    hdgs = [-10 -5 -2.5 -1.25 0 1.25 2.5 5 10]';
    % some zero values were stored as +/- eps in an older version of the gui
    data.heading(abs(data.heading)<0.01) = 0;

% remove the rest
removethese = ~ismember(data.heading,hdgs);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end

    
%% after settling on the above, run this script to generate summary data

dots3DMP_parseData


%% some plots

dots3DMP_plots


%% fit cumulative gaussians

cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
    % for probabilities, error is negative log likelihood of observing the data, which is
    % [ log(Pright(hdg)) + log(1-(~Pright(hdg))) ]
cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice))))); 

flippedGauss = @(b,hdg) 1 - ( min(max(b(1),0),1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));
    % for continuous values, error is sum squared error
flippedGauss_err = @(param,SEP,hdg) sum((flippedGauss(param,hdg)-SEP).^2);


%%

unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are unc, and plots are always search)

dots3DMP_plots_cgauss


%%
dots3DMP_plots_cgauss_forTalk






