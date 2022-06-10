function [nsEvents] = nsEventConditions(nsEvents,PDS)
% pull in PLDAPS actual conditions into Ripple Events to replace inds


%% check that unique trial numbers match completely, if not, select NS entries to match
% 04/20/2022 there should only be one NS file per PDS file!
% this will be the case if one trellis file corresponds to more than one
% PDS file.
% if so, we only want the recording trials corresponding to the desired PDS
% file

clear PDSutn
for t=1:length(PDS.data)
    PDSutn(t,:) = PDS.data{t}.unique_number;
end
NSutn = nsEvents.pldaps.unique_trial_number;
if iscell(NSutn)
    NSutn = cell2mat(NSutn');
end

matchTrials = ismember(NSutn,PDSutn,'rows');
fprintf('%d trials found in NSevents for this paradigm\n',sum(matchTrials))


% if sum(matchTrials)<length(matchTrials)
%     fprintf('%d trials did not match between PDS and NS, removing these...\n',sum(~matchTrials))
% end

% this is still a problem for RFMapping with the HeadingTheta and
% HeadingPhi! should these just be cells to make things simpler? and same
% with unique_trial_number
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

if ~isequal(PDSutn,NSutn)
    disp('unique tr nums still do not match between PDS and NS...something more serious is wrong!\n'); keyboard
end

%% 

% refactor trial data from individual cells into matrix for easier comparison with nsEvents
clear temp
% should be fixed across trials within a paradigm
fnames = fieldnames(PDS.conditions{1}.stimulus);

for t = 1:length(PDSconditions)
    
    for f=1:length(fnames)
        temp.(fnames{f})(:,t) = PDSconditions{t}.stimulus.(fnames{f});
    end
    
    % datapixx time is what is sent, use this rather than PDS.data{t}.trstart!
    PDStrstart(t) = PDSdata{t}.datapixx.unique_trial_time(1);
    
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
% across trials, but already seems super small.

% fo = fit(PDStrstart',nsEvents.Events.trStart','poly1');
% fprintf('PDS/dpixx trStart and nsEvents trstart offset is ~%.3fms',fo.p2*1e3)

% would need to go through and shift ALL the timed events
%timingEvents = {'trStart',''};
% trStartsrefit = fo(nsEvents.Events.trStart');

    
   


    

