function [nsEvents] = nsEventConditions(nsEvents,PDS,par)
% [nsEvents] = NSEVENTCONDITIONS(nsEvents,PDS)
%
% pulls in useful trial info from relevant PLDAPS file into nsEvent struct
% 
% since events are sent to Ripple via binary DigIn, conditions except
% modality are sent as indices in condition list. Here we pull the actual
% info from PLDAPS
%
% SJ updates
% 05-2022 fixed issues with mismatched trials
% 06-2022 added 'goodtrial' from PDS, since this isn't in all Ripple files
% 08-2022 added import of 1-targ behavior from PLDAPS, also not sent to Trellis

%% check that unique trial numbers match completely, if not, select NS entries to match
% note that if multiple PDS files correspond to one Trellis file, this
% will return just the trials in nsEvents which correspond to the selected
% PDS file. i.e. probably unwise to overwrite nsEvents in the wrapper!

for t=1:length(PDS.data)
    PDSutn(t,:) = PDS.data{t}.unique_number;
end
NSutn = nsEvents.pldaps.unique_trial_number;
if iscell(NSutn)
    NSutn = cell2mat(NSutn');
end

matchTrials = ismember(NSutn,PDSutn,'rows');
% fprintf('%d trials found in NSevents for this paradigm\n',sum(matchTrials))
% if sum(matchTrials)<length(matchTrials)
%     fprintf('%d trials did not match between PDS and NS, removing these...\n',sum(~matchTrials))
% end

% 06-2022
% breakfix vector was 1 entry short in one case, which messes up the
% matchTrials stuff. let's just remove it, and just add in the goodtrial logical
% later. (this needs to be fixed at some point)
try nsEvents.Events = rmfield(nsEvents.Events,'breakfix'); catch; end

% remove fields only relevant to RFmapping
if ~strcmpi(par,'RFmapping')
    try nsEvents.Events = rmfield(nsEvents.Events,'dotOrder'); catch; end
    try nsEvents.Events = rmfield(nsEvents.Events,'stimOn_all'); catch; end
    try nsEvents.Events = rmfield(nsEvents.Events,'stimOff_all'); catch; end
end


fnames = fieldnames(nsEvents.Events);
for f = 1:length(fnames)
    nsEvents.Events.(fnames{f}) = nsEvents.Events.(fnames{f})(matchTrials);
end
fnames = fieldnames(nsEvents.pldaps);
for f = 1:length(fnames)
    if strcmp(fnames{f},'unique_trial_number') && ~iscell(nsEvents.pldaps.unique_trial_number) % SJ 05-2022 some early files used a cell instead of mat, can re-create info files to fix this
        nsEvents.pldaps.unique_trial_number = nsEvents.pldaps.unique_trial_number(matchTrials,:);
    else
        nsEvents.pldaps.(fnames{f}) = nsEvents.pldaps.(fnames{f})(matchTrials);
    end
end

PDSdata = PDS.data;
PDSconditions = PDS.conditions;

NSutn=NSutn(matchTrials,:);

% in case PDS started before Trellis
matchTrials2 = ismember(PDSutn,NSutn,'rows');
PDSutn=PDSutn(matchTrials2,:);

if ~isequal(PDSutn,NSutn)
    disp('unique tr nums still do not match between PDS and NS...something more serious is wrong!\n'); keyboard
end

PDSconditions = PDSconditions(matchTrials2);
PDSdata = PDSdata(matchTrials2);

%% 
% remove fields which are all nans, irrelevent for this paradigm
% (presumably)
fnames = fieldnames(nsEvents.Events);
for f = 1:length(fnames)
    if ~iscell(nsEvents.Events.(fnames{f})) && all(isnan(nsEvents.Events.(fnames{f})))
        try nsEvents.Events = rmfield(nsEvents.Events,fnames{f}); catch; end
    end
end

% refactor trial data from individual cells into matrix for easier insertion into nsEvents
% should be fixed across trials within a paradigm
fnames = fieldnames(PDS.conditions{1}.stimulus);

clear temp
for t = 1:length(PDSconditions)
    
    for f=1:length(fnames)
        temp.(fnames{f})(:,t) = PDSconditions{t}.stimulus.(fnames{f});
    end
    
    % datapixx time is what is sent, use this rather than PDS.data{t}.trstart!
    PDStrstart(t) = PDSdata{t}.datapixx.unique_trial_time(1);
    
    nsEvents.Events.goodtrial(t) = PDSdata{t}.behavior.goodtrial;
    
    % SJ added 08-23-2022
    if strcmp(nsEvents.pldaps.parName{t},'dots3DMP')
        nsEvents.Events.oneTargChoice(t) = PDSdata{t}.behavior.oneTargChoice;
        nsEvents.Events.oneTargConf(t)   = PDSdata{t}.behavior.oneTargConf;
    end

end

% just add in the fields from PDS to nsEvents, without any further checks (except
% for the ones done above)
fnames = fieldnames(PDSconditions{1}.stimulus);
for f=1:length(fnames)
    if any(strcmp(fnames{f},{'headingTheta','hdgOrder','headingPhi'}))
        nsEvents.Events.(fnames{f}) = mat2cell(temp.(fnames{f})',ones(sum(matchTrials),1))';
    else
        nsEvents.Events.(fnames{f}) = temp.(fnames{f});
    end
end



%-----

% this is the old version, it has more checks and tests, but won't work for
% RFmapping, so needs adjustment..
%{
fnames = fieldnames(PDSconditions{1}.stimulus);
for f=1:length(fnames)
    
    % indices of the stimulus condition in nsEvents and temp should match,
    % but need to loop over blocks because list of indices may be different
    % for different blocks, even of the same paradigm
    
    % this doesn't work for RFmapping, oh and it assumes that all condition
    % are presented at least once!
    nsInds = nan(size(nsEvents.pldaps.blockNum));
    for b=1:length(nsEvents.hdr.pldaps_filenames)
        I = nsEvents.pldaps.blockNum==b;
        [~,~,nsInds(I)] = unique(temp.(fnames{f})(I));
    end
    
    if isfield(nsEvents.Events, ([fnames{f} 'Ind']))
        
        % this will also replace NaNs in Ripple Events with condition from PDS
        % (NaN in Ripple Events on brfix trials, but PDS condition is
        % specified from the beginning of the trial)
        F = nsEvents.Events.([fnames{f} 'Ind']);
        F(isnan(F)) = nsInds(isnan(F));
        
        if isequal(F,nsInds)
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
%}


% PDStrstart = PDStrstart-PDStrstart(1)+nsEvents.Events.trStart(1);

% check relative timing of PDS trial start and Ripple trial start
% could use 'linearinterp' if there is any sense that shift is different
% across trials, but already seems super small (<1ms)

% fo = fit(PDStrstart',nsEvents.Events.trStart','poly1');
% fprintf('PDS/dpixx trStart and nsEvents trstart offset is ~%.3fms',fo.p2*1e3)

% would need to go through and shift ALL the timed events
%timingEvents = {'trStart',''};
% trStartsrefit = fo(nsEvents.Events.trStart');

    
   


    

