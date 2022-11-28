% first-pass at reading in and plotting data from Marton et al.,
% "Validating a dimension of doubt in decisionmaking: A proposed
% endophenotype for obsessive-compulsive disorder"
% CF & KT 09-2021


clear

% change this to the folder where you have the data/code
cd '/Users/chris/Documents/MATLAB/Projects/Nestadt-Khushboo'

data = readtable('RDMPilotRDM_GIANT_VC_Repository.csv');

% use 'unique' on the independent variables, both to make sure we
% understand what's in the data, and to use as loop index variables below
uID = unique(data.ID) % Subject Identification #
uAGELT50 = unique(data.AGELT50) %Age less than 50 years: 1= age>50, 2=age<50
uGroupx2 = unique(data.Groupx2) % Disorder group: 1 = control, 2 = OCD 
uTest = unique(data.Test) % 0=practice, 1=$$feedback, 2=confidence, 3=$$penalty
uCoh = unique(data.Coherence) % motion strength (% coherence)

% also sanity-check dependent variables:

% first RT:
figure; hist(data.ReactionTime);
% Looks reasonable, but pretty long tail, may need to cut it off somewhere,
% ie it's possible that RTs approaching 10000 ms should not be included,
% because they may indicate the subject was not paying attention or failed 
% to respond..

% now Correct:
unique(data.Correct) % 1=incorrect, 2=correct
% Correct has some NaNs (Not a Number) which are not treated as unique by
% the function 'unique', making the above line useless. 
% These NaNs presumably indicate trials that were invalid for some reason,
% so we'll need to exclude them:
goodTrials = ~isnan(data.Correct); % (isnan generates a logical array, like '==NaN' would,
                                   % and the '~' reverses it ('NOT' operator))
% now we can try unique again:
unique(data.Correct(goodTrials)) % that's more like it

% lastly, Conf:
unique(data.Confidence) % conf rating 1-7
unique(data.Confidence(~isnan(data.Confidence))); % (showing off the flexibility of logical arrays)
% Confidence also has many NaNs, but that's because Conf judgments were
% only made in the 'test2' experiment. Let's check if that's true:
unique(data.Confidence(strcmp(data.Test,'test2')))
% aha, still some NaNs left in there. probably invalid trials. how many?
sum(isnan(data.Confidence(strcmp(data.Test,'test2'))))
% sounds about right. So we can exclude those too:
badTrials_conf = isnan(data.Confidence) & strcmp(data.Test,'test2');
goodTrials(badTrials_conf) = false;


%% now we can start looping
    % note the '%%' divides your code into 'cells' which can be run separately
    % using the 'Run Section' button

 % preallocate
pctCor = nan(length(uCoh),1);
RTmean = nan(length(uCoh),1);
RTse = nan(length(uCoh),1);
confMean = nan(length(uCoh),1);
confSE = nan(length(uCoh),1);
N = nan(length(uCoh),1);

for c = 1:length(uCoh)
    
    % start w basic plots of %corr, RT, and conf as a function of Coherence
    I = goodTrials & data.Coherence==uCoh(c);
    N(c) = sum(I);
    pctCor(c) = sum(data.Correct(I)==2) / N(c);
    
    RTmean(c) = mean(data.ReactionTime(I));
    RTse(c) = std(data.ReactionTime(I))/sqrt(N(c)); % standard approximation of standard error of the mean (SEM),
                                                      % used for plotting error bars and some statistical tests

    % I'll leave conf for you to fill in here!

    
    
end

% standard error for a percentage has a convenient formula and can be done
% outside the loop:
pctCorSE = sqrt( (pctCor.*(1-pctCor)) ./ N );


%% now make some plots!







%% next you can try to separate by participant, by age, or etc.
% For participant, follow the same logic as the coherence loop above, ie:
for n = 1:length(uID)
    
    % and you can have nested loops, like this:
    for c = 1:length(uCoh)
        
        
        
    end
    
end


% for Age, you could just add a logical statement to the 'I' index variable:
% I = goodTrials & data.Coherence==uCoh(c) & data.AGELT50==1;
% or
% I = goodTrials & data.Coherence==uCoh(c) & data.AGELT50==2;






