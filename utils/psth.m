function [Yt, x, Y, durs] = psth(spktimes,alignEvent,varargin)
% [Yt, x, Y, durs] = PSTH(spktimes,alignEvent,varargin)
%
% This code generates peri-event time histogram for each trial for a given unit, aligned
% to a single event. It also returns the total activity in the same
% interval, and duration of these interval, so can be used simultaneously,
% or separately, to calculate average firing rates in an interval.
% 
% Possible uses:
% 1. use a single vector in alignEvent to align spikes
% to a fixed duration (tEnd-tStart) around one event (e.g. motionOnset)
% e.g. [Yt, x] = psth(spktimes, [motionOn],'tStart',-1,'tEnd',2,'binSize',0.1);
%
% 2. use a Ntrx2 array in alignEvent to align spikes to a fixed time (tStart) 
% % before one event to a fixed time (tEnd) after a second event. 
% % Time '0' will be the time of the first column in alignEvent.
% In this case, the psths will have variable durations on each trial (so
% there will be NaNs in the matrix for trials shorter than the max dur).
% e.g. [Yt, x] = psth(spktimes, [motionOn saccOnset],'tStart',-1,'tEnd',1,'binSize',0.1);
%
% 3. use a Ntrx2 array in alignEvent to calculate average firing rate
% between two events
% [~,~,Y,durs] = psth(spktimes, [motionOn saccOnset],'tStart',0,'tEnd',0);
% 
% 
% the code should work regardless of time units (s, ms etc), as long as
% they are consistent for all variables! Defaults for tStart, tEnd, and
% binSize are all in seconds.
%
%
% INPUTS: 
% required arguments
%   spktimes   - vector of spike times [nspk x 1]
%   alignEvent - times of task events for aligning spikes to
%         [ntr x 1, or ntr x 2] 
%   
% optional arguments
%       tStart, tEnd  - time relative to alignEvent to start and end psth
% if alignEvent is Ntrx1, tmin and tmax will both be relative to that single alignEvent
% if alignEvent is Ntrx2, then tStart will be relative to first alignEvent, and tEnd will be relative to second
%       binSize     - the size of each spike count bin, default 50ms
%       alltrials   - if 0, only include trials after the first spike and before the final spike (inclusive)
%       normalize   - if 1, divide by binsize to give firing rate over time, and
%       divide overall 'y' count by durs to get average firing rate
%
% OUTPUTS:
%       Yt   - firing rate or spike counts over time [Ntrs x nbins]
%       x    - time of Yt bins (center of bin) relative to alignment event [nbins x 1]
%       Y    - mean firing rate in interval [Ntrs x 1]
%       durs - duration of intervals [Ntrs x 1]
%
% 
% SJ 08-2022 created
% SJ 09-2022 added functionality for variable duration psths, and aligning
% to first column in alignEvent, regardless of order of events
% 
% TODO: 
%
% use some flags to skip binning sections if user only wants average rate
% in interval

%% parse the inputs

p = inputParser;
p.addRequired('spktimes');
p.addRequired('alignEvent');

var_names = {'tStart','tEnd','binSize','alltrials','normalize'};
defaults = {-1,1,0.05,false,true};

for v=1:length(var_names)
    p.addParameter(var_names{v},defaults{v});
end

parse(p,spktimes,alignEvent,varargin{:});
p = p.Results;

%% deal with the alignEvent var

Ntr=length(alignEvent);
if (Ntr==0), disp('No good trials...'); return; end

% TODO: allow user to specify explicitly which column to align to?
if size(alignEvent,2)==2
    [alignEvent,EvOrder] = sort(alignEvent,2);
    whichEventAlign = EvOrder(1); % align to whichever event was originally first column in alignEvent
else     % if alignEvent is a vector, tStart and tEnd are both relative to same event              
    alignEvent(:,2) = alignEvent(:,1);
    whichEventAlign = 1;
end

% get the start and end times from alignment event matrix
tr_starts = alignEvent(:,1) + p.tStart;
tr_ends   = alignEvent(:,2) + p.tEnd;

durs   = tr_ends - tr_starts;

%% get indices of first and last trials that we actually want
itr_start = 1; itr_end   = Ntr;

if ~p.alltrials && ~isempty(spktimes)
    [~,t] = min(abs(tr_starts-spktimes(1)));
    if ~isempty(t), itr_start=t; end
    [~,t] = min(abs(tr_ends-spktimes(end)));
    if ~isempty(t), itr_end=t; end
end

%%  define the x-axis 'edges' for PSTH, and pre-allocate fr

% to get the overall extent of bins, we need to re-calculate the start OR
% end relative to the desired alignment event, since tStart and tEnd are
% initially set relative to separate events in alignEvent.
% (if alignEvent has just one column, this is redundant, but still fine)

if whichEventAlign==2
    tStarts_new = tr_starts - alignEvent(:,2);
    tStart_new  = min(tStarts_new);
    tEnd_new    = p.tEnd;
elseif whichEventAlign==1
    tEnds_new   = tr_ends - alignEvent(:,1);
    tEnd_new    = max(tEnds_new);
    tStart_new  = p.tStart;
end

if p.tStart<0

    % this way, the event you are aligning to has a bin on either side, it
    % does not fall within a bin 
    x = 0:-p.binSize:tStart_new; % backwards from zero
    x1= 0:p.binSize:tEnd_new; % forwards from zero

    % deal with the ends
    if abs(x(end)-tStart_new)>0.01, x(end+1)=x(end)-p.binSize; end
    if abs(x1(end)-tEnd_new)>0.01, x1(end+1)=x1(end)+p.binSize; end

    % for some reason, tEnd_new-x1 has some small non-zero component
    % if x(end)~=p.tStart, x(end+1)=x(end)-p.binSize; end
    % if x1(end)~=tEnd_new, x1(end+1)=x1(end)+p.binSize; end

    % stitch them together now
    x=[fliplr(x) x1(2:end)];

else
    % if both tStart_new and tEnd_new are >0, then just do it the simple way
    x = tStart_new:p.binSize:tEnd_new;
end


nBins = length(x)-1;
Yt = nan(Ntr,nBins);
Y  = nan(Ntr,1);

% now compute the psth within dur limits on valid trials
for itr=itr_start:itr_end
   if isnan(alignEvent(itr,whichEventAlign)), continue; end
   
   spkInd = spktimes >= tr_starts(itr) & spktimes <= tr_ends(itr);
   
   Y(itr) = sum(spkInd);    
   if Y(itr)==0, Yt(itr,:) = 0; continue; end
   
   index1=spktimes(spkInd)-alignEvent(itr,whichEventAlign);
   
   Yt(itr,:) = histcounts(index1,x);

   % if we have two events (with variable delay) in alignEvent, discard the
   % bins on trials which extend beyond the dur on that trial 
   % if aligning to earlier event,  bins will be at the back end of the trials. 
   % if aligning to the later event, bins will be at the start.
   if whichEventAlign==1
       [~,endPos] = min(abs(x-tEnds_new(itr)));
       Yt(itr,endPos:end) = NaN;
   elseif whichEventAlign==2
       [~,startPos] = min(abs(x-tStarts_new(itr)));
       Yt(itr,1:startPos) = NaN;
   end
   
end

% shift 'x' from edges of bins to middle of bins, helpful for later plotting
x=x+diff(x(1:2))/2;
x(end)=[]; % cut off the end, so that x is nBins long


if p.normalize  
    Yt = Yt / p.binSize; % normalize spike count by bin width to get firing rate!
    Y = Y ./ durs;       % normalize total spike count by total interval time
end