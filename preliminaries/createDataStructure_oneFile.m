function data = createDataStructure_oneFile(PDS,subject)

%% createDataStructure_oneFile(PDS)
% create data struct using just one PDS struct
% mostly for testing purposes
% 08-20-21 SJ

data = struct;
data.choice = []; % initialize this one field, you'll see why

fieldExcludes = {'leftEarly','tooSlow','fixFP','FPHeld','eyeXYs','corrLoopActive','goodtrial', ...
    'timeTargDisappears','probOfMemorySaccade','leftTargR','leftTargTheta', ...
    'rightTargR','rightTargTheta','audioFeedback','textFeedback','rewardDelay'};


T = length(data.choice); % set trial counter to continue where left off (0, to start)
for t = 1:length(PDS.data) % loop over trials for this file,
    if isfield(PDS.data{t}.behavior,'choice') % and save out the data, excluding trials with
        % missing data (choice is a good marker for this)
        T = T+1; % increment trial counter
        
        % independent variables are stored in PDS.conditions.stimulus
        fnames = fieldnames(PDS.conditions{t}.stimulus);
        fnames(ismember(fnames,fieldExcludes)) = [];
        for F = 1:length(fnames)
            eval(['data.' fnames{F} '(T,1) = PDS.conditions{t}.stimulus.' fnames{F} ';']);
        end
        % EXCEPT! dot duration, and anything else in
        % PDS.data.stimulus, i.e. things that are generated
        % on each trial and not pre-configured. For now
        % will need to hard-code them here.
        if isfield(PDS.data{t}.stimulus,'dotDuration')
            data.duration(T,1) = PDS.data{t}.stimulus.dotDuration;
        end
        if isfield(PDS.data{t}.stimulus,'dotPos')
            data.dotPos{T,1} = PDS.data{t}.stimulus.dotPos;
        end
        
       
        
        % dependent variables are stored in PDS.data.behavior
        fnames = fieldnames(PDS.data{t}.behavior);
        fnames(ismember(fnames,fieldExcludes)) = [];
        for F = 1:length(fnames)
            % SJ 07-20, correct defaults to logical but
            % throws an error for NaN - change to
            % double
            if strcmp(fnames{F},'correct'), data.correct(T,1) = 0; end
            eval(['data.' fnames{F} '(T,1) = PDS.data{t}.behavior.' fnames{F} ';']);
        end
        
        % noticed a couple extra things we need, not in either place -CF 02-2021
        try
            data.oneTargPDW(T,1) = PDS.data{t}.postTarget.markOneConf;
        catch
        end
        try
            data.delayToPDW(T,1) = PDS.data{t}.postTarget.delayToConfidence;
        catch
        end
        
    end
end
clear PDS


%% some bookkeeping: manually merge any vars whose names have changed, etc

if strcmp(subject(1:5),'human')
    data.conf = data.saccEndPoint;
else
    % do we still need this? SJ 07/2020
    if isfield(data,'postDecisionConfidence')
        if isfield(data,'PDW')
            if length(data.postDecisionConfidence)<length(data.PDW) && length(data.PDW)==length(data.choice)
                data.PDW(1:length(data.postDecisionConfidence)) = data.postDecisionConfidence;
                data = rmfield(data,'postDecisionConfidence');
            else
                error('unsure, diagnose by looking at data');
            end
        elseif length(data.postDecisionConfidence)==length(data.choice)
            data.PDW = data.postDecisionConfidence;
            data = rmfield(data,'postDecisionConfidence');
        else
            error('unsure, diagnose by looking at data');
        end
    end
end
if isfield(data,'saccEndPoint')
    data = rmfield(data,'saccEndPoint'); % either way, this gets removed
end

if isfield(data,'oneTargTrial')
    if isfield(data,'oneTargChoice')
        if length(data.oneTargTrial)<length(data.oneTargChoice) && length(data.oneTargChoice)==length(data.choice)
            data.oneTargChoice(1:length(data.oneTargTrial)) = data.oneTargTrial;
            data = rmfield(data,'oneTargTrial');
        else
            error('unsure, diagnose by looking at data');
        end
    elseif length(data.oneTargTrial)==length(data.choice)
        data.oneTargChoice = data.oneTargTrial;
        data = rmfield(data,'oneTargTrial');
    else
        error('unsure, diagnose by looking at data');
    end
end
    
    
if isfield(data,'oneConfTargTrial')
    if isfield(data,'oneTargConf')
        if length(data.oneConfTargTrial)<length(data.oneTargConf) && length(data.oneTargConf)==length(data.choice)
            data.oneTargConf(1:length(data.oneConfTargTrial)) = data.oneConfTargTrial;
            data = rmfield(data,'oneConfTargTrial');
        else
            error('unsure, diagnose by looking at data');
        end
    elseif length(data.oneConfTargTrial)==length(data.choice)
        data.oneTargConf = data.oneConfTargTrial;
        data = rmfield(data,'oneConfTargTrial');
    else
        error('unsure, diagnose by looking at data');
    end
end
       

disp('done.');

