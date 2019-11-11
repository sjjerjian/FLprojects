%% this is the beginning of retroactive PLDAPS (here, PDS) cleanup.

% get all the files and erase all the garbage we don't use. this will save
% disk space and also make finding things easier.

% as of 08/2019, goal is no longer to modify files on the server, only
% local copies. thus we can be more brutal with our cuts. Later, if we
% really feel like it or need to free up space on the NAS, we can run this
% over there.

% these may now be some important toggles for some analyses:
removeAnalogData = 1;
removeTimingData = 1;

        % consider deleting:
% PDS.data{t}.datapixx.adc.data: what are the six channels? add an index var
    % PDS.data{t}.datapixx.eyeXind etc. has just 1-3, what are 4-6?
    
    
        % definitely delete:
% PDS.initialParameters
% PDS.initialParameterNames
% PDS.initialParametersMerged
% PDS.functionHandles
% PDS.conditionNames (redundant)
% PDS.data{t}.pldaps.lastBgColor
% PDS.data{t}.behavior.reward: unused
% PDS.data{t}.mouse: only save samples if useMouseAsEyePos option selected
% PDS.data{t}.pldaps.draw: just framerate, not interested
% PDS.data{t}.pldaps.trialStates
% PDS.data{t}.pldaps.goodTrial: move to .behavior, then remove .pldaps entirely
% PDS.data{t}.plexon: obvi
% PDS.data{t}.task.*: move to .behavior, and delete .task
% PDS.data{t}.flagNextTrial
% PDS.data{t}.iFrame
% PDS.data{t}.state
% PDS.data{t}.keyboard: maybe in the future keep it?
% PDS.data{t}.currentFrameState
% PDS.data{t}.eyeX
% PDS.data{t}.eyeY
% PDS.data{t}.pupil
% PDS.data{t}.remainingFrameTime
% PDS.data{t}.framePreLastDrawIdleCount
% PDS.data{t}.framePostLastDrawIdleCount

        % misc
% PDS.data{t}.stimulus.eyeXYs: move to .behavior
% PDS.data{t}.stimulus.dotX/Y/Z/Size: make a single NaN if vestib only
% PDS.data{t}.postTarget: check, deleted for now
% PDS.data{t}.reward : now just fixRewarded, add things



% temp, for home:
% load('/Users/chris/Documents/MATLAB/PLDAPS_data/human/human20190625dots3DMP1413.PDS','-mat');


% subject = 'humanCXD';


% first thing though is to rename a few variables:
allFiles = dir(localDir);
for d = 1:length(dateRange)
    for f = 3:length(allFiles) % skip 1+2, they are are "." and ".."
        if contains(allFiles(f).name, subject) ... % check if the target file is in localDir
           && contains(allFiles(f).name, num2str(dateRange(d))) ...
           && contains(allFiles(f).name, paradigm)
       
            disp(['loading ' allFiles(f).name]);
            load([localDir allFiles(f).name],'-mat'); % load it. this is the time limiting step;
                                                      % will eventually change how data are saved to make this faster
                               
            try PDS = rmfield(PDS,'initialParameters'); catch; end
            try PDS = rmfield(PDS,'initialParameterNames'); catch; end
            try PDS = rmfield(PDS,'initialParametersMerged'); catch; end
            try PDS = rmfield(PDS,'functionHandles'); catch; end
            try PDS = rmfield(PDS,'conditionNames'); catch; end
            
            for t = 1:length(PDS.data) % loop over trials for this file
                
                try PDS.data{t}.behavior = rmfield(PDS.data{t}.behavior,'reward'); catch; end
                if removeAnalogData
                    try PDS.data{t}.datapixx = rmfield(PDS.data{t}.datapixx,'adc'); catch; end
                end
                try PDS.data{t} = rmfield(PDS.data{t},'mouse'); catch; end
                try PDS.data{t}.behavior.goodtrial = PDS.data{t}.pldaps.goodtrial; catch;  end
                try PDS.data{t} = rmfield(PDS.data{t},'pldaps'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'plexon'); catch; end
                try PDS.data{t}.behavior.eyeXYs = PDS.data{t}.stimulus.eyeXYs; catch; end
                try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'eyeXYs'); catch; end
                try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameStimOn'); catch; end
                try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameStimOff'); catch; end
                try PDS.data{t}.behavior.timeTargDisappears = PDS.data{t}.task.timeTargDisappears; catch; end
                try PDS.data{t}.behavior.probOfMemorySaccade = PDS.data{t}.task.probOfMemorySaccade; catch; end
                try PDS.data{t} = rmfield(PDS.data{t},'task'); catch; end
                try PDS.data{t} = rmfield(PDS.data{t},'postTarget'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'flagNextTrial'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'iFrame'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'state'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'keyboard'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'currentFrameState'); catch; end                               
                if removeTimingData
                    try PDS.data{t} = rmfield(PDS.data{t},'timing'); catch; end
                end
                try PDS.data{t} = rmfield(PDS.data{t},'eyeX'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'eyeY'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'pupil'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'remainingFrameTime'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'framePreLastDrawIdleCount'); catch; end                               
                try PDS.data{t} = rmfield(PDS.data{t},'framePostLastDrawIdleCount'); catch; end                               

                % rename fields (dots3DMP)
                try PDS.conditions{t}.stimulus = renameStructField(PDS.conditions{t}.stimulus,'conflictAngle','delta'); catch; end
                try PDS.conditions{t}.stimulus = renameStructField(PDS.conditions{t}.stimulus,'stimCondition','modality'); catch; end
                try
                    if PDS.conditions{t}.stimulus.modality==1
                        PDS.data{t}.stimulus.dotX = NaN;
                        PDS.data{t}.stimulus.dotY = NaN;
                        PDS.data{t}.stimulus.dotZ = NaN;
                        PDS.data{t}.stimulus.dotSize = NaN;
                    end
                catch    
                end
                
            end
                        
            save([localDir allFiles(f).name],'PDS');
            
        end
    end
end


