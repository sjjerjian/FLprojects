function dots3DMP_NeuralStruct_runCleanUp(dataStruct,parSelect,minRate,minTrs)

fprintf('Cleaning up data, resaving\n')

dataStruct_clean = dataStruct;
removeEntireSession = false(size(dataStruct));

for s = 1:length(dataStruct)
    clear parUnits numSpikes numTrials enoughTrials

    % sessions that don't have all pars in parSelect get marked for
    % removal, but do it at the end so as not to mess up the loop counter
    if ~all(isfield(dataStruct(s).data,parSelect))
        removeEntireSession(s) = true;
        continue
    end

    for par = 1:length(parSelect)

        units  = dataStruct(s).data.(parSelect{par}).units;
        events = dataStruct(s).data.(parSelect{par}).events;

        parUnits(par,:)  = units.cluster_id;
        numSpikes(par,:) = cellfun(@length,units.spiketimes);
        numTrials(par,:) = length(events.trStart);


        stimCondList = [events.heading; events.modality; events.coherence; events.delta]';

        for u = 1:length(units.cluster_id)

            if ~isempty(units.spiketimes{u})
                [~,t] = min(abs(events.trStart-units.spiketimes{u}(1)));
                if ~isempty(t), itr_start=t; end
                [~,t] = min(abs(events.trStart-units.spiketimes{u}(end)));
                if ~isempty(t), itr_end=t; end
            end

            [uStimConds,~,ic]    = unique(stimCondList(itr_start:itr_end,:),'rows');
            [nTrConds,~]         = hist(ic,unique(ic));
            enoughTrials(par,u)  = all(nTrConds>=minTrs);
        end

    end

    numTrials       = repmat(numTrials,1,size(parUnits,2));

    parSpikeRate = numSpikes ./ numTrials;

    removeThese = any(parSpikeRate<minRate | ~enoughTrials,1);

    for par = 1:length(parSelect)
        
        units  = dataStruct(s).data.(parSelect{par}).units;
        units.cluster_id(removeThese) = [];
        units.cluster_type(removeThese) = [];
        units.spiketimes(removeThese) = [];

        % overwrite
        dataStruct_clean(s).data.(parSelect{par}).units = units;
    end
end

dataStruct_clean(removeEntireSession) = [];
dataStruct = dataStruct_clean;

file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '_neuralData_clean.mat'];
save([localDir(1:length(localDir)-length(subject)-7) file], 'dataStruct');