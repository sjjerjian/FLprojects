% generate .mat file with recording info

clear info
createBinFile = 1; % set to 0 to skip (just create nsEvents)

addpath(genpath('C:\Program Files (x86)\Ripple\Trellis\Tools'));
addpath(genpath('C:\Users\fetschlab\Documents\MATLAB\dots3DMP'));

%% START ACTUAL OFFLINE PROCESSING HERE
% one cell per recording date (penetration)

%% 2022-06-15, pen 27

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220615;
info.pen         = 27;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[3, 5]};
info.chanlist            = [3]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1502 1512 1528 1606 1626];
info.trellis_filenums    = [1 2 3 4 5];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 1 1 2 2];
info.comments            = {'lots of brfix initially, coh 1 only','1 only, 20 per cond. initially lots of low bets','1+2+3, 0.7 coh, stopped early','L bias (ves)',''};
info.chanInterest        = {[3],[3],[3],[3],[3]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9100 9100 9100 9120 9120]}; % MDI depth
info.cellcomments        = {{'very clean','consistent, probably a small MU too','','cell got a bit smaller but still decent',''}};
% this is effectively becoming chanInterest at this stage, can probably

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-05-26, pen 26

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220526;
info.pen         = 26;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[1, 0]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1255 1310 1444 1452 1528 1541];
info.trellis_filenums    = [1 1 1 1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning','RFMapping'};
info.rec_group           = [1 1 1 1 1 1];
info.comments            = {'1+2 only','all 3 mods. some R bias on high coh. Forgot to set numReps so stopped early','stuck in zero coh corr loop, need to fix this','all high bet on vis...','',''};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9400 9400 9400 9400 9400 9400]}; % MDI depth
info.cellcomments        = {{'','','','','',''}};
% this is effectively becoming chanInterest at this stage, can probably
% skip it unless there is something relevant about signal in general

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2022-05-23, pen 25

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220523;
info.pen         = 25;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[1, 1]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1214 1228 1249 1303 1325 1352 1502];
info.trellis_filenums    = [1 1 1 1 1 1 1];
info.par                 = {'RFMapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 1 1 1 1 1 1];
info.comments            = {'just 4 thetas','all mods','all mods, ves PDW drifting low','ves only','vis + comb, more high bets','tried to encourage more low bets, but overall similar','all 3 mods'};
info.chanInterest        = {[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10200,10200,10200,10200,10200,10200,10200]}; % MDI depth
info.cellcomments        = {{'','','','','','',''}};
% this is effectively becoming chanInterest at this stage, can probably
% skip it unless there is something relevant about signal in general

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

% 2022-05-20, pen 24

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220520;
info.pen         = 24;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[1, 1]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1305 1316 1348 1411 1605];
info.trellis_filenums    = [1 2 2 2 2];
info.par                 = {'RFMapping','RFMapping','dots3DMPtuning','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 2 2 2 2];
info.comments            = {'ignore, testing mode','2 apertures, 2R, 8 theta','all 3 mods','all 3. more low bets on ves, especially at the start',''};
info.chanInterest        = {[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9900 9600 9600 9600 9600]}; % MDI depth
info.cellcomments        = {{'','','','',''}};
% this is effectively becoming chanInterest at this stage, can probably
% skip it unless there is something relevant about signal in general

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-17, pen 23

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220517;
info.pen         = 23;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[2, 5]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1151 1158 1206 1237 1249 1335 1353 1402];
info.trellis_filenums    = [1 1 1 1 1 1 1 1];
info.par                 = {'VesMapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMPtuning','VesMapping'};
info.rec_group           = [1 1 1 1 1 1 1 1];
info.comments            = {'no FP, horiz plane only','1+2,0.7 coh only','1 only','2+3, mostly betting high','changed rew/pens','lots more low bets, esp on vis','just 0.7 coh again','no FP'};
info.chanInterest        = {[],[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9100 9100 9100 9100 9100 9100 9100 9100]}; % MDI depth
info.cellcomments        = {{'ch1,ch30 clear units, + others','ch31, weirdly regular firing','huge cell on ch10','ch7 got smaller, large sharp unit on ch21 + 32','2 cells on 5, huge on 16+24, 7 got smaller','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-13, pen 22

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220513;
info.pen         = 22;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1152 1207 1220 1251 1350 1358 1406];
info.trellis_filenums    = [1 1 1 1 1 1 1];
info.par                 = {'VesMapping','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMP','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 1 1 1 1];
info.comments            = {'no fix req, but FP shown.','1+2','0.05 1-targ, Ves only. Good','Vis+Comb, stopped early as increasingly betting high','Still almost all PDW high','Mostly high bet','More low bets, changed incentives'};
info.chanInterest        = {[],[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[8900 8900 8900 8900 8900 8900 8900]}; % MDI depth
info.cellcomments        = {{'huge cell on 23/24, stayed throughout','','','','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-12, pen 21

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220512;
info.pen         = 21;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[3, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1238 1252 1525];
info.trellis_filenums    = [1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','RFMapping'};
info.rec_group           = [1 1 1];
info.comments            = {'','All 3 mods. brfix a lot, R biased on choices, PDW good','15 trials, brfix a lot, stopped early'};
info.chanInterest        = {[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[8100 8100 8100]}; % MDI depth
info.cellcomments        = {{'','lots of activity! Ref switch was the issue...',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-05-05, pen 20

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220505;
info.pen         = 20;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'10034_32ch_65um_80mm'};
info.gridxy              = {[2, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1441 1454 1517];
info.trellis_filenums    = [1 1 1];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1];
info.comments            = {'no IMU rec. 0.2+0.7coh','2 only','1+3, PDW steeper',};
info.chanInterest        = {[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[10100 10100 10100]}; % MDI depth
info.cellcomments        = {{'lots of active chs but not clear SU it seems','activity decreased a bit? ch11 seems to be task-responsive',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-08, pen 19

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220408;
info.pen         = 19;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[6, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1345 1420 1435 1506 1645];
info.trellis_filenums    = [1 2 3 4 5];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 2 2 2 2];
info.comments            = {'','','','','',''};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[6300 6500 6500 6500 6500]}; % MDI depth
info.cellcomments        = {{'','','','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-07, pen 18

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220407;
info.pen         = 18;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[6, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1305 1325 1338 1402 1414 1445];
info.trellis_filenums    = [1 2 3 4 5 6];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 2 2 2];
info.comments            = {'1+2 only','1+2 only','Vis only. Betting high, aborted','1+2 only','1 only, good!','2+3, good! 0 hdgs always rewarded right, BUG!!, 0.1 high target hidden'};
info.chanInterest        = {[],[],[],[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[7200 7200 7200 7700 7700 7700]}; % MDI depth
info.cellcomments        = {{'few units on mapped chs 12,21,27. Probably going for reward','','12 got smaller','12-31 activity','v clear on 12 and 21, 27 gets smaller','12,21,22,31 all good'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-05, pen 17

clear info
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220405;
info.pen         = 17;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[2, 3]};
info.chanlist            = [1:32]; % include all chs here, remove dead chs at kilosort chanMap stage

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1311 1325 1421];
info.trellis_filenums    = [1 2 3];
info.par                 = {'dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1];
info.comments            = {{'1+2 only, 0.2+0.7 coh','Vis only,0.2+0.7,5% 1-target PDW, betting low ~10%. Slight R bias','Ves + Comb, 0.2+0.7. ~10% low bet again, Low and High Coh Comb look similar'}};
info.chanInterest        = {[],[],[]};        % chans of interest for online thresholded spikes

% >=1 value, different across probes (recording info)
info.depths              = {[9400 9400 9400]}; % MDI depth
info.cellcomments        = {{'ch30 possibly more than one unit','ch5 as well, ch30 maybe bigger until the end','ch30 smaller again, definitely task responsive'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,createBinFile);

%% 2022-04-01, pen 16
% 1 value, fixed across session
clear info

info.subject     = 'lucio';
info.date        = 20220401;
info.pen         = 16;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'DBC D-A 32ch'};
info.probe_ID            = {'8714_32ch_65um_80mm'};
info.gridxy              = {[4, 2]};
info.chanlist             = [1:32];

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1406];
info.trellis_filenums    = [1];
info.par                 = {'dots3DMPtuning'};
info.rec_group           = [1];
info.comments            = {{'1+2 only, 0.7 coh only'}};
info.chanInterest        = {[]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[7005]} ; % MDI depth
info.cellcomments        = {{'nothing really there'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-14, pen 14
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220324;
info.pen         = 14;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[3, 0]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1432,1456];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 2];
info.comments            = {{'all mods, 0.7 coh only','all mods, 0.7 coh only'}};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[9929,10632]}; % MDI depth
info.cellcomments        = {{'very active, multiple cells. probably MU','also very active, probably MU'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0); % no bin file for single ch recording

%% 2022-03-14, pen 11
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220314;
info.pen         = 11;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridxy              = {[2,-2]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1400,1528];
info.trellis_filenums    = [1 2];
info.par                 = {'RFMapping','dots3DMPtuning'};
info.rec_group           = [1 2];
info.comments            = {{'tried to sample space, 30+20apertureD','all mods, ~80 trials then forgot to stop Trellis recording!'}};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[4966,7924]}; % MDI depth
info.cellcomments        = {{'very sharp, clear cell eventually, clear response to stim','probably MU,wider'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-11, pen 10 
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220311;
info.pen         = 10;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[6,3]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1315,1522,1540,1551];
info.trellis_filenums    = [1 2 3 4];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 2 2 2];
info.comments            = {{'4trs per condition, 0.2 low coh','full set, started second PLDAPs at end, ignore last few trials (combined)','maybe ignore','started recording after PLDAPS'}};
info.chanInterest        = {[3],[3],[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths for different recordings

% >=1 value, different across probes (recording info)
info.depths              = {[6076,11733,11733,11870]}; % MDI depth
info.cellcomments        = {{'good cell','','','cell better isolated, same as before?'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-08, pen 9
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220308;
info.pen         = 9;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[4,4]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1441,1502];
info.trellis_filenums    = [1,2];
info.par                 = {'dots3DMPtuning','dots3DMP'};
info.rec_group           = [1 1];
info.comments            = {{'full protocol','short block'}};
info.chanInterest        = {[3],[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[12454,12454]}; % MDI depth
info.cellcomments        = {{'very clean cell','initially good, cell response changed a bit later'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-07, pen 8
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220307;
info.pen         = 8;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'FHC'};
info.gridxy              = {[2,6]};
info.chanlist             = 3;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1424];
info.trellis_filenums    = [1];
info.par                 = {'dots3DMPtuning'};
info.rec_group           = [1];
info.comments            = {{'90 trials'}};
info.chanInterest        = {[3]};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[8591]}; % MDI depth
info.cellcomments        = {{'spike sorted online, ok. one mid-sized unit, and a larger one, possibly responding to reward'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-03-01, pen 5
% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220301;
info.pen         = 5;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type          = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridxy              = {[2,3]};
info.chanlist             = 1;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1425, 1445, 1505, 1519, 1540, 1614];
info.trellis_filenums    = [1 2 3 4 5 7];
info.par                 = {'dots3DMPtuning','dots3DMPtuning','dots3DMPtuning','dots3DMP','dots3DMP','dots3DMPtuning'};
info.rec_group           = [1 2 3 3 3 4];
info.comments            = {{'','','','extra +1.5 hdg,vis+comb, mostly betting high','ves only, good',''}};  
info.chanInterest        = {1,1,1,1,1,1};        % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)  
info.depths              = {[5099, 5496, 6079, 6079, 6079, 7506]}; % MDI depth
info.cellcomments        = {{'not much, some cell bursts','one big cell low FR, got smaller, another smaller one','','cell very stable','got a bit smaller','huge cell'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);

%% 2022-02-23, pen 4

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220223;
info.pen         = 4;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type           = {'Single electrode'};
info.probe_ID            = {'AlphaOmega'};
info.gridtype             = 'odd_sym';
info.gridxy               = {[1,1]};
info.chanlist             = 1;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1438, 1451, NaN, 1515, NaN, 1600];
info.trellis_filenums    = 1:6;
info.par                 = {'dots3DMPtuning','RFMapping','RFMapping','RFMapping','dots3DMPtuning','dots3DMPtuning'};
info.rec_group           = [1 1 2 3 4 5];
info.comments            = {{'','','ignore','forgot to stop the Trellis recording for a while','ignore',''}}; 
info.chanInterest        = {1,1,1,1,1,1};         % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[7583, 7583, NaN, 7910, NaN, 8719]}; 
info.cellcomments        = {{'probably lost SU, then got another one?','','lost the cell','','',''}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
save([info.filepath savefilename],'info','-mat')
processTrellisData(info,1);

%% 2022-02-18, pen 3

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220218;
info.pen         = 3;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS
% info.filepath = ['C:\Users\fetschlab\Trellis\dataFiles\' num2str(info.date) '\'];

% >=1 value, fixed across session
info.probe_type           = {'Single electrode'};
info.probe_ID             = {'FHC'};
info.gridxy               = {[0,0]};
info.chanlist             = 1;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [1622,NaN];
info.trellis_filenums    = [1 2];
info.par                 = {'dots3DMP','RFMapping'};
info.rec_group           = [1 2];
info.comments            = {{'R biased choices','pldaps recording wasn''t set up'}};
info.chanInterest        = {1,1};         % chans of interest for online thresholded spikes, can be different lengths

% >=1 value, different across probes (recording info)
info.depths              = {[9772,9521]}; % MDI depth
info.cellcomments        = {{'maybe an SU, bad online sorting, not super responsive','clearer unit'}};

savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);


%% TESTING , FAKE RECORDINGS
if 0

% 1 value, fixed across session
info.subject     = 'lucio';
info.date        = 20220615;
info.pen         = NaN;
info.gridtype    = 'odd_sym';

info.filepath = ['\\172.30.3.33\homes\fetschlab\data\' info.subject '\' info.subject '_neuro\' num2str(info.date) '\']; % save directly to NAS

% >=1 value, fixed across session
info.probe_type          = {'NaN'};
info.probe_ID            = {'NaN'};
info.gridxy              = {[NaN,NaN]};
info.chanlist             = 1:32;

% >=1 value, fixed across probes (task info)
info.pldaps_filetimes    = [NaN];
info.trellis_filenums    = [1 1 1 2];
info.par                 = {'RFMapping','dots3DMPtuning','dots3DMP','dots3DMP'};
info.rec_group           = [1 1 1 2];
info.comments            = {{''}};
info.chanInterest        = {[],[]};        % chans of interest for online thresholded spikes, can be different lengths for different recordings

% >=1 value, different across probes (recording info)
info.depths              = {[NaN]}; % MDI depth
info.cellcomments        = {{'',''}};

% savefilename = sprintf('%s%ddots3DMP_info.mat',info.subject,info.date);
% save([info.filepath savefilename],'info','-mat')
processTrellisData(info,0);
end