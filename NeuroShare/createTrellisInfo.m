% generate .mat file with recording info
% at the end of each info cell, call processTrellisData
% this will 
%   a) create a RippleEvents structure containing timings of all
%       key task events, trial conditions, and behavioral outcomes
%   b) add a .spkData field to this events structure containing any
%   threshold crossings and spikes detected online (through Trellis hoops)
%   c) save raw continuous data (.ns5 format) as int16 binary for offline
%   sorting in Kilosort (one file per probe/electrode)
% creation of nsEvents struct requires createEventStruct_dots3DMP_NS, and
% nested getData_NS


% based on info.probe_type and probe_ID


addpath(genpath('C:\Program Files (x86)\Ripple\Trellis\Tools'));
addpath(genpath('C:\Users\fetschlab\Documents\MATLAB\dots3DMP'));

%% QUICK PROCESS (minimum version for semi-online use e.g. tuningCurves)

if 0
    info.subject     = 'test';
    info.date        = 20220315;
    info.probe_type          = {'test'};
    info.probe_ID            = {'test'};
    info.chanlist             = {1};

    info.pldaps_filetimes    = [1232 1250];
    info.trellis_filenums    = [1 2];
    info.par                 = {'dots3DMPtuning','VesMapping'};
    info.rec_group           = [1 2];
    info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

    % >=1 value, different across probes (recording info)
    info.depths              = {[NaN NaN]}; % MDI depth
    info.chanInterest        = {{1,1}};        % chans of interest for online thresholded spikes, can be different lengths
    processTrellisData(info);
end


%% START ACTUAL OFFLINE PROCESSING HERE
% one cell per recording date (penetration)

%% 2022-03-14 
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220314;
info.pen         = 11;
info.gridtype    = 'odd_sym';

info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridxy              = {[2,-2]};
info.chanlist             = {3};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1400,1528];
info.trellis_filenums    = [1 2];
info.par                 = {'RFMapping','dots3DMPtuning'};
info.rec_group           = [1 2];

% >=1 value, different across probes (recording info)
info.depths              = {[4966,7924]}; % MDI depth
info.chanInterest        = {{3,3}};        % chans of interest for online thresholded spikes, can be different lengths
info.comments            = {{'tried to sample space, 30+20apertureD','all mods, ~80 trials then forgot to stop Trellis recording!'}};
info.cellcomments        = {{'very sharp, clear cell eventually, clear response to stim','probably MU,wider'}};

savefilename = sprintf('%s%d dots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info);

%% 2022-03-11 
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220311;
info.pen         = 10;
info.gridtype    = 'odd_sym';

info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[6,3]};
info.chanlist             = {3};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1315,1522,1540,1551];
info.trellis_filenums    = [1 2 3 4];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 2 2 2];

% >=1 value, different across probes (recording info)
info.depths              = {[6076,11733,11733,11870]}; % MDI depth
info.chanInterest        = {{3,3,3,3}};        % chans of interest for online thresholded spikes, can be different lengths
info.comments            = {{'4trs per condition, 0.2 low coh','full set, started second PLDAPs at end, ignore last few trials (combined)','maybe ignore','started recording after PLDAPS'}};
info.cellcomments        = {{'good cell','','','cell better isolated, same as before?'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info);

%% 2022-03-08 
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220308;
info.pen         = 9;
info.gridtype    = 'odd_sym';

info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[4,4]};
info.chanlist             = {3};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1441,1502];
info.trellis_filenums    = [1,2];
info.par                 = {'dots3DMPtuning','dots3DMP'};
info.rec_group           = [1 1];

% >=1 value, different across probes (recording info)
info.depths              = {[12454,12454]}; % MDI depth
info.chanInterest        = {{3,3}};        % chans of interest for online thresholded spikes, can be different lengths
info.comments            = {{'full protocol','short block'}};
info.cellcomments        = {{'very clean cell','initially good, cell response changed a bit later'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info);

%% 2022-03-07 
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220307;
info.pen         = 8;
info.gridtype    = 'odd_sym';

info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[2,6]};
info.chanlist             = {3};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1424];
info.trellis_filenums    = [1];
info.par                 = {'dots3DMPtuning'};
info.rec_group           = [1];

% >=1 value, different across probes (recording info)
info.depths              = {[8591]}; % MDI depth
info.chanInterest        = {{3}};        % chans of interest for online thresholded spikes, can be different lengths
info.comments            = {{'90 trials'}};
info.cellcomments        = {{'spike sorted online, ok. one mid-sized unit, and a larger one, possibly responding to reward'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info);

%% 2022-03-01 
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220301;
info.pen         = 5;
info.gridtype    = 'odd_sym';

info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridxy              = {[2,3]};
info.chanlist             = {1,1};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1425, 1445, 1505, 1519, 1540, 1614];
info.trellis_filenums    = [1 2 3 4 5 7];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 2 3 3 3 4];

% >=1 value, different across probes (recording info)
info.depths              = {[5099, 5496, 6079, 6079, 6079, 7506]}; % MDI depth
info.chanInterest        = {{1,1,1,1,1,1}};        % chans of interest for online thresholded spikes, can be different lengths
info.comments            = {{'','','','extra +1.5 hdg,vis+comb=, mostly betting high','ves only, good',''}};  
info.cellcomments        = {{'not much, some cell bursts','one big cell low FR, got smaller, another smaller one','','cell very stable','got a bit smaller','huge cell'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info);

%% 2022-02-23

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220223;
info.pen         = 4;
info.gridtype    = 'odd_sym';

info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type           = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridtype             = 'odd_sym';
info.gridxy               = {[1,1]};
info.chanlist             = {1,1};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1438, 1451, NaN, 1515, NaN, 1600];
info.trellis_filenums    = 1:6;
info.par                 = {'dots3DMPtuning','RFMapping','RFMapping','RFMapping','dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 1 2 3 4 5];

% >=1 value, different across probes (recording info)
info.depths              = {[7583, 7583, NaN, 7910, NaN, 8719]}; 
info.chanInterest        = {{1,1,1,1,1,1}};         % chans of interest for online thresholded spikes, can be different lengths
info.comments            = {{'','','ignore','forgot to stop the Trellis recording for a while','ignore',''}}; 
info.cellcomments        = {{'probably lost SU, then got another one?','','lost the cell','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info);

%% 2022-02-18

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220218;
info.pen         = 3;
info.gridtype    = 'odd_sym';
info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type           = {'Single electrode'};
info.probe_ID             = {'FHC'};
info.gridxy               = {[0,0]};
info.chanlist             = {1,1};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1622,NaN];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMP','RFMapping'};
info.rec_group           = [1 2];

% >=1 value, different across probes (recording info)
info.depths              = {[9772,9521]}; % MDI depth
info.chanInterest        = {{1,1}};         % chans of interest for online thresholded spikes, can be different lengths
info.comments            = {{'R biased choices','pldaps recording wasn''t set up'}};
info.cellcomments        = {{'maybe an SU, bad online sorting, not super responsive','clearer unit'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info);

if 0
%% TESTING , FAKE RECORDINGS

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220306;
info.pen         = NaN;
info.gridtype    = 'odd_sym';

info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'NaN'};
info.probe_ID            = {'NaN'};
info.gridxy              = {[NaN,NaN]};
info.chanlist             = {NaN};

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1528];
info.trellis_filenums    = [3];
info.par                 = {'VesMapping'};
info.rec_group           = [1];
% >=1 value, different across probes (recording info)
info.depths              = {[NaN NaN]}; % MDI depth
info.chanInterest        = {{NaN NaN}};        % chans of interest for online thresholded spikes, can be different lengths for different recordings
info.comments            = {{'',''}};
info.cellcomments        = {{'',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info);
end