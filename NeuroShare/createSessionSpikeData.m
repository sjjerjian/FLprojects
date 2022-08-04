

%% ENTER INFO HERE

subject = 'lucio';
date = 20220311;
spikeSorter = 'mkSort'; % 'kilosort'

%%
dataFolder = 'C:\Users\fetschlab\Trellis\dataFiles';
thisSession = fullfile(dataFolder,num2str(date));
cd(thisSession)

load([subject num2str(date) 'dots3DMP_info.mat'])

for f=1:length(info.trellis_filenums)

    file = sprintf('%s%ddots3DMP%04d',subject,date,info.trellis_filenums(f));

    % load in the events file, we'll add this in alongside the spikes data
    load([file '_RippleEvents.mat']) 

    % bit of clean-up
    events = nsEvents.Events;
    fnames = fieldnames(nsEvents.pldaps);
    for fn=1:length(fnames)
        events.(fnames{fn}) = nsEvents.pldaps.(fnames{fn});
    end

    cd(sortfolder)
        
    switch spikeSorter
        case 'mkSort'
            spkfiles = dir('waveforms_*.mat');

            for s = 1:length(spkfiles)
                load(spkfiles(s).name);

                spikeTimes = waveforms.spikeTimes;
                spikeWaveforms = waveforms.alignedWaves;

                unitIDs = unique(waveforms.units);
                sortquality = waveforms.ratings.ratings;

                clear spks
                for u=1:length(unitIDs)
                    iu = unitIDs(u);
                    if iu>0
                        spks{iu}.unitnum = iu;
                        spks{iu}.quality = sortquality(iu);

                        inds = waveforms.units==iu;

                        spks{iu}.spikeTimes = spikeTimes(inds);
                        spks{iu}.spikeWaves = waveforms.alignedWaves(inds);
                        spks{iu}.nSpikes    = sum(inds);
                    end
                end
            end
        case 'kilosort'

    end

    cd .. 

    % do we want to save each individual recording as a separate file, even
    % if it's the same cell across paradigms, or different settings? need
    % to consider the options here
%     save([file '_spks.mat'],'spks','events','info');




% todo
% kilosort implementation
% consider layout for adding different 'paradigms' together, or if e.g.
% different modalities of 3DMP are recorded separately for the same
% units...

end

