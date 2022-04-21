%{
SJ 04/2022
Notes on the neural data processing pipeline for the Ripple Grapevine/MP rig.


%============
1. RECORDING AND RAW DATA

Neural data and timestamps of key behavioral events are recorded on Trellis
system. 

Timestamps and online sorted digital spike events are stored in .nev files
1kHz analog data (eye signals and IMU accelerometer signals) are in .ns2
30kHz neural data (raw broadband, lfp, and high-pass filtered) are in .ns5

These raw data files should be copied to NAS.
%============
2. PRE-PROCESSING

1) input recording notes into a cell block within createTrellisInfo.m for
the given recording day (see template). This will create and save an info
struct containing all the info for the day's session.

2) processTrellisData will then be run using this info file, to create an
nsEvents struct (saved in [subject][date]dots3DMP[filenum]_RippleEvents.mat), and a binary int16 file for
Kilosort (saved as [subject][date]_[recording set number]. 
There will be one binary file per recording 'set' (1 or more Ripple Events), so that
waveforms from the same recording depth across different PLDAPS experiments
can be sorted together, saved as _1, _2 etc. See *SPIKE SORTING* note below

nsEvents
- contains trial condition and outcome information, and timing of key
behavioral events (in seconds)

since trial conditions are sent as binary codes, all except modality are
coded as the index in condition list. offline recovery of actual condition
value is done using nsEventConditions.m during post-processing.

- for online window discriminated spikes (if any), nsEvents will contain a spkData field with
spike timing and cluster information

- for offline alignment purposes, nsEvents will contain a field called
analogData, with important header information about the neural data. For
aligning spikes to behavioral events, spike times extracted from Kilosort need to be corrected by the
'startTime' in the relevant nsEvents.analogData for aligning to events in
nsEvents, or the nsEvents can be corrected by the timestamps in order to
show data from multiple recordings together. See *SPIKE SORTING* note below

%============
3. NEXT STEP PROCESSING

1. Offline spike sorting for probe recordings should be done with Kilosort,
pointing to the relevant binary file for a given date and recording set.
(main_kilosort.m, SJ parameters on local Trellis computer, with appropriate config files set up)
After manual curation, necessary .phy output files will be in the _[set]
folder where the binary file is located (on NAS).

For single electrode recordings, .ns5 files can be loaded into MKsort in
Matlab, and sorted manually.

2. creating neural data structure
for a given subject, dateRange, and paradigm, run
dots3DMP_NeuralPreProcessing, which calls createNeuralDataStruct to
generate a structure containing one row per recorded unit.


*** NOTE ON SPIKE SORTING ***

Generally it is better for Kilosort to handle recordings done at the same
depth (with presumably the same neurons) as one file, instead of a user
trying to manually reconcile clusters across separately sorted files.

On the other hand, it seems more straightforward to make sure each Trellis
recording is linked to just one PLDAPS file on the experimental computer.
(This means stopping and starting Trellis recordings for each paradigm e.g.
tuning and task, even at the same depth).
(Why is it more straightforward?)
- 1. Because we need to reconcile the PLDAPS conditions list with the
nsEvents one (and also nexonar potentially)
- 2. The PLDAPS events are different for different paradigms, so
constructing nsEvents is more parsimonious / cleaner
- 3. We can decide offline how to combine recordings for
analysis/visualization e.g. two dots3DMP recordings for different modalities
- 4. 

The alternative is to just run the Trellis recording for the whole time we
are recording different 'paradigms' at the same location, and deal with the
events offline. So we create one binary file for a given recording set for
Kilosort purposes, and store the 'start' and 'end' time of each recording
within a set for offline re-alignment.
