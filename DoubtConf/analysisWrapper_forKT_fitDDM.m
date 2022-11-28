% first-pass at reading in and plotting data from Marton et al.,
% "Validating a dimension of doubt in decision-making: A proposed
% endophenotype for obsessive-compulsive disorder"
% CF & KT 09-2021

clear
close all

% change this to the folder where you have the data/code
% cd '/MATLAB Drive'
cd '/Users/chris/Documents/MATLAB/Projects/DoubtConf(Nestadt)'
data = readtable('RDMPilotRDM_GIANT_VC_Repository.xlsx');
data = table2struct(data,'ToScalar',true);

% use 'unique' on the independent variables, both to make sure we
% understand what's in the data, and to use as loop index variables below
uID = unique(data.ID); % Subject Identification #
uAGELT50 = unique(data.AGELT50); %Age less than 50 years: 1= age>50, 2=age<50
uGroupx2 = unique(data.Groupx2); % Disorder group: 1 = control, 2 = OCD 
uTest = unique(data.Test);       % 0=practice, 1=$$feedback, 2=confidence, 3=$$penalty
uCoh = unique(data.Coherence);   % motion strength (% coherence)

goodTrials = ~isnan(data.Correct);
badTrials_conf = isnan(data.Confidence) & strcmp(data.Test,'test2');
goodTrials(badTrials_conf) = false;


    
%% prep data to fit DDM

% select subset of data with all three measures, and convert to format
% expected by fitting code:
data_orig = data; clear data;
ind = goodTrials & ~isnan(data_orig.Confidence) & ~isnan(data_orig.ReactionTime) & data_orig.ReactionTime<=9000; % exclude a small number of very long RT trials

data.ageLT50 = data_orig.AGELT50(ind);
data.group = data_orig.Groupx2(ind);

data.correct = data_orig.Correct(ind)-1; % was 1=err, 2=corr, change to 0:1
data.coherence = data_orig.Coherence(ind);

% because motion direction is not indicated in the data, assign a random
% direction to each trial
dirs = [0 180]; % 0=left, 180=right
data.direction = randsample(dirs,sum(ind),'true');
data.direction = data.direction';

% the fitting code requires coherence to be 'signed', with negative values
% indicating leftward and positive indicating rightward
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);

% choice is 1 for rightward and 0 for leftward, which we can assign based
% on data.correct and the direction we randomly assigned above
data.choice = nan(size(data.correct));
data.choice(data.correct==1 & data.direction==0) = 1;
data.choice(data.correct==1 & data.direction==180) = 0;
data.choice(data.correct==0 & data.direction==0) = 0;
data.choice(data.correct==0 & data.direction==180) = 1;

data.RT = data_orig.ReactionTime(ind)/1000; % convert to seconds
data.conf = data_orig.Confidence(ind);
data.dur = 9000 * ones(size(data.conf));

% fitting functions currently all require binary confidence (ie PDW), so
% we need to median-split the rating scale data, for now
data.PDW = data.conf >= median(data.conf);

save doubtconf.mat
disp('done');
    