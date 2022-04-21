function [nsEvents] = nsEventConditions(nsEvents,PDS)

% pull in PLDAPS actual conditions into Ripple Events to replace inds
% could just pass in the conditions lists, but double use this to clean-up any issues in the
% data, check alignment of timings and potentially add more info in Ripple
% time)

% before this function, will need checks to determine which PLDAPS file and
% which trellis file to use (i.e. based on info file)

% clear;clc
% cd /Users/stevenjerjian/Desktop/FetschLab/Analysis/data/20220311
% load('lucio20220311dots3DMPtuning1315.mat')
% load('lucio20220311dots3DMP0001_RippleEvents.mat')



% check that unique trial numbers match completely, if not, select NS entries to match
% i.e. remove any NS trials that occurred before PDS (this is only possible
% if PLDAPS was started twice within one recording)
% 04/20/2022 there should only be one NS file per PDS file!
for t=1:length(PDS.data)
    PDSutn(t,:) = PDS.data{t}.unique_number;
end
NSutn = cell2mat(nsEvents.pldaps.unique_trial_number');

matchTrials = all(ismember(NSutn,PDSutn),2);

if sum(matchTrials)<length(matchTrials)
    fprintf('%d trials did not match between PDS and NS, removing these...',sum(~matchTrials))
end

fnames = fieldnames(nsEvents.Events);
for f = 1:length(fnames)
    nsEvents.Events.(fnames{f}) = nsEvents.Events.(fnames{f})(matchTrials);
end
fnames = fieldnames(nsEvents.pldaps);
for f = 1:length(fnames)
    nsEvents.pldaps.(fnames{f}) = nsEvents.pldaps.(fnames{f})(matchTrials);
end

% get conditions and data for Ripple Events trials (should be 1:nTrials anyway, unless
% Ripple recording missed a trial e.g. PLDAPS started first)
PDSdata = PDS.data(nsEvents.pldaps.iTrial);
PDSconditions = PDS.conditions(nsEvents.pldaps.iTrial);

% re-compute these now
clear PDSutn
for t=1:length(PDSdata)
    PDSutn(t,:) = PDSdata{t}.unique_number;
end
NSutn = cell2mat(nsEvents.pldaps.unique_trial_number');

if length(PDSconditions)~=length(nsEvents.Events.modality)
    disp('number of trials does not match!\n'); keyboard
elseif ~isequal(PDSutn,NSutn)
    disp('unique tr nums still do not match between PDS and NS...something more serious is wrong!\n'); keyboard
end


fnames = fieldnames(PDSconditions{1}.stimulus);

for t = 1:length(PDSconditions)

    for f=1:length(fnames)
        temp.(fnames{f})(t) = PDSconditions{t}.stimulus.(fnames{f});
    end
    
    % do something with data too? get RTs etc
%     fnames = fieldnames(PDS.data{1}.behavior);

%     for f=1:length(fnames)
%         temp2.(fnames{f})(t) = PDSdata{t}.behavior.(fnames{f});
%     end

    % datapixx time is what is actually sent to Ripple, don't use
    % PDS.data{t}.trstart!
    PDStrstart(t) = PDSdata{t}.datapixx.unique_trial_time(1);  
    
    
end



for f=1:length(fnames)
    
    % indices of the stimulus condition in nsEvents and temp should match
    [~,~,ic] = unique(temp.(fnames{f}));
    
    if isfield(nsEvents.Events, ([fnames{f} 'Ind']))
        
        % replace NaNs in Ripple Events with condition from PDS
        F = nsEvents.Events.([fnames{f} 'Ind']);
        F(isnan(F)) = ic(isnan(F));
        
        if isequal(F,ic')
            nsEvents.Events.(fnames{f}) = temp.(fnames{f});
        else
            fprintf('%s indices do not match between PDS and NS\n',fnames{f})
            keyboard
        end
        
    elseif strcmp(fnames{f},'modality')
        
        F = nsEvents.Events.modality;
        nsEvents.Events.modality(isnan(F)) = temp.modality(isnan(F));
        
        if ~isequal(nsEvents.Events.modality,temp.modality)
            fprintf('modality does not match between PDS and NS\n')
            keyboard
        end
        
    end
end

% PDStrstart = PDStrstart-PDStrstart(1)+nsEvents.Events.trStart(1);

% check relative timing of PDS trial start and Ripple trial start
% could use 'linearinterp' if there is any sense that shift is different
% across trials, but already seems super small.

% fo = fit(PDStrstart',nsEvents.Events.trStart','poly1');
% fprintf('PDS/dpixx trStart and nsEvents trstart offset is ~%.3fms',fo.p2*1e3)

% would need to go through and shift ALL the timed events
%timingEvents = {'trStart',''};
% trStartsrefit = fo(nsEvents.Events.trStart');

    
   


    

