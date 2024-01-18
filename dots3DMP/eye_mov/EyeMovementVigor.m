% Confidence-conditioned Eye Movement Vigor Analysis
% for Yueh-Chen
% Oct 2023


% 1. Visualize eye movements for 'good' trials for 4 trial outcomes (or 2
% outcomes - collapse across right and left trials)

% 2. define threshold/windows in X (or X+Y) to calculate velocity of
% left/right saccade

% 3. compare velocity for high vs low bets
% hypothesis: choice saccade velocity for high bets should be higher than for low bets

% 4. compare velocity for correct vs error trials (across all bets, and/or
% within high and low bets separately)

% 5. see if you can identify any change of mind trials (trajectory going to
% one target then turning round - either in choice or wager). 
% saccades are ballistic so this is likely
% to be very rare, may only see a few in the larger dataset)

close all
clc

%% define variables
% extract a fixed time window around the time when eyes leave fixation point
% this way, we can create a matrix where every trial is a fixed length
tmin = -0.2; % seconds
tmax = 0.6;
sample_rate = 0.001;
% the number of individual you want to plot
N = 200;
% set thresholds for velocity
r_thr = 0.3;
l_thr = 0.2;
%how to charactrize the data sets
right_left = 1;
high_low = 1;
correct_error = 0;


%% keep good trials only
goodtrial = ~isnan(data.choice) & ~isnan(data.PDW) & data.oneTargConf==0;

% extract data
[em_data_mat, em_time_vec] = extractgoodtrial(tmin,tmax,sample_rate,data,goodtrial);


%% sort the trial
% think about right left/ high low/ correct and error

% find right/left, and high/low
[r_hi, r_lo, l_hi, l_lo] = sorttrial(data, goodtrial,em_data_mat);

%% plot individual trials and mean

% plot mean
plotmean(em_data_mat,em_time_vec,r_hi, r_lo,l_hi,l_lo)
%plot individual
plotindividualtrials(N, r_hi,r_lo,l_hi, l_lo,em_data_mat,em_time_vec);

%% Velocity

[v_r_hi] = velocity(r_hi, em_data_mat, sample_rate,r_thr, N);
[v_r_lo] = velocity(r_lo, em_data_mat, sample_rate,r_thr, N);
[v_l_hi] = velocity(l_hi, em_data_mat, sample_rate,l_thr, N);
[v_l_lo] = velocity(l_lo, em_data_mat, sample_rate,l_thr, N);

%% t-test and ploting

[p_r, p_l] = ttestandplot(v_r_hi, v_r_lo, v_l_hi,v_l_lo);
fprintf("ttest results:\nRight: %f\nLeft: %f\n", p_r, p_l);

%% functions

function [em_data_mat, em_time_vec] = extractgoodtrial(tmin,tmax,sample_rate,data,goodtrial)

em_time_vec = tmin:sample_rate:tmax;

% extract variables from each trial
rt_time = data.timeGoCue(goodtrial); % time that eyes leave fixation point (for creating matrix)
em_data = data.adc.data(goodtrial);  % datapixx data (chs 1 and 2 are eye X and Y position)
em_time = data.adc.time(goodtrial);  % datapixx time (for x-axis)

% pre-allocate size of big matrix
em_data_mat = zeros(3, length(em_time_vec), length(em_time));

% ...and fill it with the data
for i = 1:length(em_time)
    [~, pos(1)] = min(abs(em_time{i}-rt_time(i)-tmin));
    [~, pos(2)] = min(abs(em_time{i}-rt_time(i)-tmax));

    temp = em_data{i}(1:3, pos(1):pos(2));

    if all(abs(temp(1,:)) < 0.8) && ...
       all(abs(temp(1,1:(-tmin/sample_rate-40))) < 0.2) && ...
       all(abs(temp(1,(-tmin/sample_rate+40):end)) > 0.4)
  
    em_data_mat(:, 1:size(temp, 2), i) = temp; % store trial data in big matrix
    continue
    end
    em_data_mat(1, 1, i) = NaN;
end


% subtract starting value so all traces start from zero
em_data_mat = em_data_mat - mean(em_data_mat(:, 1:100, :), 2);

% moving average, to smooth data a bit
% use with care, because an acausal averaging might obfuscate some timing
% effects which are key to calculating threshold crossings)
% em_data_mat = movmean(em_data_mat, 10, 2);

% NOTE: eye tracker flips Y axis, flip it back (so that it looks like the screen)
em_data_mat(2, :, :) = -em_data_mat(2, :, :);
end

function [r_hi, r_lo, l_hi, l_lo] = sorttrial(data, goodtrial,em_data_mat)
% & data.heading==0 
% to select for heading
r_hi_good = find(data.PDW(goodtrial)==1 & data.choice(goodtrial)==2);
r_lo_good = find(data.PDW(goodtrial)==0 & data.choice(goodtrial)==2);
l_hi_good = find(data.PDW(goodtrial)==1 & data.choice(goodtrial)==1);
l_lo_good = find(data.PDW(goodtrial)==0 & data.choice(goodtrial)==1);

r_hi_idx = ~isnan(em_data_mat(1,1,r_hi_good)); 
r_lo_idx = ~isnan(em_data_mat(1,1,r_lo_good)); 
l_hi_idx = ~isnan(em_data_mat(1,1,l_hi_good)); 
l_lo_idx = ~isnan(em_data_mat(1,1,l_lo_good));

r_hi = r_hi_good(r_hi_idx);
r_lo = r_lo_good(r_lo_idx);
l_hi = l_hi_good(l_hi_idx);
l_lo = l_lo_good(l_lo_idx);

fprintf("trial number:\nright high %d\nright low %d\nleft high %d\nrleft low %d\n",length(r_hi), length(r_lo), length(l_hi),length(l_lo));
end

function [V] = velocity(idx,em_data_mat, sample_rate, thr, N)

len = length(idx);
V = zeros(len,1);
P = sort(randperm(len,N));
j = 1;
figure
for i = 1:len

    pos_x = em_data_mat(1,:,idx(i));
    pos = find(abs(pos_x)>thr, 1, 'first' );

    %set threshold for position
%     hi_thr = mean(abs(pos_x(pos+20:pos+120)))-std(abs(pos_x(pos+20:pos+120)));
%     lo_thr = mean(abs(pos_x(pos-120:pos-20)))+std(abs(pos_x(pos-120:pos-20)));
% 
%     start = find(abs(pos_x(1:end))<(lo_thr), 1, 'last' );
%     stop = find(abs(pos_x(1:end))>hi_thr, 1, 'first' );

    % set threshold for velocity
    vec_x = abs(movmean(diff(pos_x),10));
    vec_thr = mean(vec_x(1:100))+2*std(vec_x(1:100));
    start = find(abs(vec_x(1:pos))<(vec_thr), 1, 'last' );
    stop = pos+ find(abs(vec_x(pos:end))<(vec_thr), 1, 'first' );
    
    if ~isempty(start) && ~isempty(stop)
       
        V(i,1) = abs((pos_x(stop)-pos_x(start))/ ((stop-start+1)*sample_rate));
        
        if j <= N
        if i == P(j)
            
            subplot(2,2,2)
            pos_x_stop = abs(pos_x(stop-40:stop+10))-mean(abs(pos_x(stop:stop+20)));
            plot(pos_x_stop);
            hold on
    
            subplot(2,2,1)
            plot(pos_x(start-10:start+40));
            hold on

            subplot(2, 2, 3);
            plot(vec_x(start-9:start+41));
            hold on

            subplot(2, 2, 4);
            vec_x = abs(movmean(diff(pos_x),10));
            plot(vec_x(stop-39:stop+11));
            hold on

            j = j+1;
        end
        end
    end
end

subplot(2,2,1)
title("start, P(x)");
set(gca,'xlim', [1 51]);
plot([10,10], ylim, 'r--');
subplot(2,2,2)
title('stop, P(x)');
set(gca,'xlim', [1 51]);
plot([40,40], ylim, 'r--');
subplot(2,2,3)
title('start, Vx(t)');
set(gca,'xlim', [1 51]);
plot([10,10], ylim, 'r--');
subplot(2,2,4)
title('stop, Vx(t)');
set(gca,'xlim', [1 51]);
plot([40,40], ylim, 'r--');
end

function [p_r, p_l] = ttestandplot(v_r_hi, v_r_lo, v_l_hi,v_l_lo)

maxSampleSize = max([numel(v_r_hi), numel(v_r_lo), numel(v_l_hi), numel(v_l_lo)]);

v_r_hi = [v_r_hi; NaN(maxSampleSize - numel(v_r_hi), 1)];
v_r_lo = [v_r_lo; NaN(maxSampleSize - numel(v_r_lo), 1)];
v_l_hi = [v_l_hi; NaN(maxSampleSize - numel(v_l_hi), 1)];
v_l_lo = [v_l_lo; NaN(maxSampleSize - numel(v_l_lo), 1)];

data = [v_r_hi, v_r_lo, v_l_hi, v_l_lo];
labels = {'Right High', 'Right Low', 'Left High', 'Left Low'};

figure;
boxplot(data, 'Labels', labels, 'Widths', 0.5);
title('Eye Movement Vigor');
ylabel('Velocity');

[~,p_r] = ttest2(v_r_hi, v_r_lo);
[~,p_l] = ttest2(v_l_hi, v_l_lo);

end

function plotmean(em_data_mat,em_time_vec,r_hi, r_lo,l_hi,l_lo)
figure;

yl = [1, 0.5, 0.5];
titles = {'X', 'Y'};

for ch = 1:2
    subplot(2, 2, ch); hold on
    plot(em_time_vec, squeeze(mean(em_data_mat(ch, :, r_hi), 3)), 'r', 'linew', 1)
    plot(em_time_vec, squeeze(mean(em_data_mat(ch, :, r_lo), 3)), 'm', 'linew', 1)
    
    plot(em_time_vec, squeeze(mean(em_data_mat(ch, :, l_hi), 3)), 'b', 'linew', 1)
    plot(em_time_vec, squeeze(mean(em_data_mat(ch, :, l_lo), 3)), 'c', 'linew', 1)

    set(gca,'ylim', [-1 1]*yl(ch))
    title(titles{ch})
end

subplot(2, 2, [3, 4]); hold on; title('X-Y')
plot(squeeze(mean(em_data_mat(1, :, r_hi), 3)), squeeze(mean(em_data_mat(2, :, r_hi), 3)), 'r', 'linew', 1)
plot(squeeze(mean(em_data_mat(1, :, r_lo), 3)), squeeze(mean(em_data_mat(2, :, r_lo), 3)), 'm', 'linew', 1)
plot(squeeze(mean(em_data_mat(1, :, l_hi), 3)), squeeze(mean(em_data_mat(2, :, l_hi), 3)), 'b', 'linew', 1)
plot(squeeze(mean(em_data_mat(1, :, l_lo), 3)), squeeze(mean(em_data_mat(2, :, l_lo), 3)), 'c', 'linew', 1)
axis([-1 1 -1 1]); axis square;
end

function plotindividualtrials(N,r_hi,r_lo,l_hi, l_lo,em_data_mat,em_time_vec)

figure;

% plot random subset of trials
r_hi = r_hi(randperm(length(r_hi), N));
r_lo = r_lo(randperm(length(r_lo), N));
l_hi = l_hi(randperm(length(l_hi), N));
l_lo = l_lo(randperm(length(l_lo), N));


yl = [1, 0.5, 0.5];
titles = {'X', 'Y'};

for ch = 1:2
    subplot(2, 2, ch); hold on
    plot(em_time_vec, squeeze(em_data_mat(ch, :, r_hi)), 'r', 'linew', 1)
    plot(em_time_vec, squeeze(em_data_mat(ch, :, r_lo)), 'm', 'linew', 1)
    
    plot(em_time_vec, squeeze(em_data_mat(ch, :, l_hi)), 'b', 'linew', 1)
    plot(em_time_vec, squeeze(em_data_mat(ch, :, l_lo)), 'c', 'linew', 1)

    set(gca,'ylim', [-1 1]*yl(ch))
    title(titles{ch})
end

subplot(2, 2, [3, 4]); hold on; title('X-Y')
plot(squeeze(em_data_mat(1, :, r_hi)), squeeze(em_data_mat(2, :, r_hi)), 'r', 'linew', 1)
plot(squeeze(em_data_mat(1, :, r_lo)), squeeze(em_data_mat(2, :, r_lo)), 'm', 'linew', 1)
plot(squeeze(em_data_mat(1, :, l_hi)), squeeze(em_data_mat(2, :, l_hi)), 'b', 'linew', 1)
plot(squeeze(em_data_mat(1, :, l_lo)), squeeze(em_data_mat(2, :, l_lo)), 'c', 'linew', 1)
axis([-1 1 -1 1]); axis square;
end




