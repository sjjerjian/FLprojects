%% SJ 2021-07-19

%% load the data

clear; clc; close all
load('sampleNexonarData_forArielle.mat')

%% Useful information

%{
 This file contains data storing the recorded movement of the motion platform via a Nexonar
 Motion Visualizer on individual trials of the cue-combination paradigm.
 There are two structs in the file; hdr and nex
%}

disp(hdr) % hdr (header) is for data management purposes only, it contains some filename
            % information from the processing pipeline.
%
disp(nex) % THIS IS THE USEFUL ONE
 
%{
 The original data is stored in one .txt file for the entire recording.
 Here, this .txt file has already been parsed to extract the motion data and
 information for each trial in the recording.

 nex contains four further structs:
    -1-  .nexdata
        
        Each cell of .nexdata contains the data from one trial e.g.
        nex.nexdata{1} is the nexonar data from trial 1

        Columns 1 and 2 are relative and absolute timestamps
        Columns 3-5 are [x,y,z] positions of the platform, in absolute
        co-ordinates
        (Note that 'heading direction' in the context of the task is
        specified by the x, i.e. column 3)

    -2-  .pldaps
        iTrial (trial number), and trialSeed for each trial, to match with
        experimental computer. Can be ignored for now.

    -3-  .conditions
        contains the stimulus conditions for each trial (stored as 1xn
        vectors, where each entry corresponds to matching entry in nexdata)
        MODALITY 
            1,2,3 (vestibular, visual, combined)
        HEADING 
             the heading angle (-ve to the left, +ve to the right)
        COHERENCE 
            proportion of coherently moving dots in the optic flow stimulus. 
            for vestibular trials, coherence is set by default to the low one, although this doesn't mean anything
        
        DELTA 
            cue conflict angle in combined condition. Zero means no conflict, 
            +ve means visual to the right, vestibular to the left (relative to heading angle in HEADING)
            -ve means vestibular to the right, visual to the left. 
            (e.g. heading angle +3 and delta -2 means vestibular is at +4 and visual is at +2)
            Some datasets don't have any conflict trials. All single-cue conditions default to delta = 0.

    -4-  .behavior
        GOODTRIAL (not the same as correct)
            1 = successfully completed trial = 1, 0 = incomplete/aborted trial (break fixation)
            if goodtrial = 0, other fields (choice, correct, RT, conf) should be NaN
            for most analyses, remove/ignore these trials
        CHOICE 
            1 = left, 2 = right
        CORRECT
            1 = correct, 2 = error 
        RT
            time of response relative to start of motion, in seconds
        conf
            post-decision wager confidence report
            1 = high confidence, 2 = low confidence

%}

%%

% ***** TO DO *****
% (optional) - generate psycometric curves to show the proportion of
% rightward choices as a function of heading angle (nex.conditions.heading), for different
% stimulus conditions (modality, coherence, delta)
%
% - create plots of motion trajectories for good trials, and split by
% heading/modality (individual trials and averages)
% assess variance in motion trajectories over time for fixed stimulus
% consider splitting trials by RT (short vs long, or multiple bins)
%
% plot proportion of righward choices using actual motion instead of
% 'heading' condition (you would likely need to compute some average
% measure over time e.g. deviation from zero), and bin data into intervals
% to have enough trials in each interval. Not sure if this will work in
% practice, but should be able to compare it to the normal psychometric
% curves!
%
% focusing on heading == 0, what are the deviations?
% if there is deviation away from zero on some/all trials, does this in
% anyway relate to the left/right choices the monkey makes?

% are there any deviations in motion trajectories for given headings which
% influence confidence judgements or RT?

% Note that some motion trajectories may look weird
% some of these are due to breakfix (aborted trials)
% others may be because issues around the response time
% breakfix trials (goodtrial == 0) should be ignored from analyses

% you may need to subtract a baseline from the traces (so that the first 1
% or n timepoints average is zero)

% ***** Additional notes *****
%
% currently I don't have the timestamps from the experimental computer
% stored here so it isn't possible to indicate when response times occur on
% individual trials (the timestamps in nex.nexdata do not line up with the RT variable) - I am working to add this in

% - if you notice anything weird in the data, don't hesitate to ask/mention
% it...the recording and pre-processing is still a work in progress so any
% comments are helpful!

% similarly if you think of something you would like to show or do, but don't have the
% necessary info (like the RT timestamps point above), please mention it - we can also add things to the
% recording/saving pipeline to have this information available

