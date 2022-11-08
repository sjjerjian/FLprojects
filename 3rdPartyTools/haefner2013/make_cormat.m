function cormat = make_cormat 

% population_likelihood.m 
% Trial-by-trial simulation of the behavioral experiment, computing
% likelihood functions by drawing spike counts from real neuronal data.

% THIS VERSION SIMULATES RESPONSES WITH A DESIRED NOISE CORRELATION LEVEL,
% ASSUMING GAUSSIAN NOISE WITH A VARIANCE = FANO*MEAN

close all
clear all
% BASE_PATH = 'C:\Program Files\MATLAB\R2007a\work\';
% BASE_PATH = '/home/chris/work/'; % if running on linux cluster
BASE_PATH = '/Users/crfetsch/Documents/DeAngelis_Angelaki/work/'; % mac

% batch_name = 'all_cells_yos'
% batch_name = 'all_cells_web'
% batch_name = 'all_cells'

% batch_name = 'cong_yos_gt0'
% batch_name = 'cong_web_gt0'
% batch_name = 'cong_both_gt0'

% batch_name = 'cong_yos_gt0.4'
% batch_name = 'cong_web_gt0.4'
batch_name = 'cong_both_gt0.4'

% batch_name = 'rocable';
% batch_name = 'opp_both_0.4';

% batch_name = 'failed_F_test';
% batch_name = 'failed_AIC';

% batch_name = 'Fig7_black_bars';
% batch_name = 'weights_really_dont_change';

load([BASE_PATH 'condition_list_neuro.mat']);
% load([BASE_PATH 'condition_list_interp.mat']);

if BASE_PATH(1)=='/'
    load([BASE_PATH 'neurons_indiv_trials/' batch_name '.mat']);
    files = textread([BASE_PATH 'neurons (temp)/' batch_name '.txt'],'%q');
else
    load([BASE_PATH 'neurons_indiv_trials\' batch_name '.mat']);
    files = textread([BASE_PATH 'neurons (temp)\' batch_name '.txt'],'%q');
end

num_reps = 200;
num_sim_repeats = 1;
save_data = 1;

use_logL = 0;
use_fake_data = 1;
    make_fake = 0;
    
% noise_corr_max = 0;
noise_corr_max = 0.1503; % trained
% noise_corr_max = 0.2683; % naive

noisemodel = 'gausian'; fano = 1.5;
% noisemodel = 'poisson'; fano = 1;

repetition = num_reps;
cell_num = length(resp_mean_all);
% method = 'spline';
method = 'linear';

wichmann = 0;
plotflag = 1; % for psychometrics
prior = 0;
Jan = 0;
Yong = 0;
FILE = ['SIM'];

stimtype = condition_list(:,1);
coherence = condition_list(:,2);
heading = condition_list(:,3);
delta = condition_list(:,4);

unique_stim_type = unique(stimtype);
unique_coherence = unique(coherence);
unique_heading = unique(heading);
unique_delta = unique(delta);
unique_amplitude = .30;
unique_azimuth = 90;
con_txt = 'delta';
con_txt_2 = 'coherence';

x = -10:0.1:10;
total_tr = length(stimtype) * num_reps;
stim = zeros(1,total_tr); coh = zeros(1,total_tr);
hdg = zeros(1,total_tr); delt= zeros(1,total_tr);
choice = zeros(1,total_tr);

% sort cells by 'pref dir' (slope around straight ahead)
for n = 1:length(resp_mean_all)
    Y = resp_mean_all{n}{3,2,2}';
    X = [ones(length(unique_heading),1) unique_heading];
    b = regress(Y,X);
    slope(n) = b(2);
end
[y,sortind] = sort(slope);
resp_mean_all = resp_mean_all(sortind);

%------------------------------------------------
% construct covariance matrix, separately for each stim type.
% method from Shadlen et al. 1996;  code modified from Yong's.

% first convert mean firing rates to Yong's format, and *interpolate*
% X = unique_heading';
X = [-10 -3.5 -1.3 0 1.3 3.5 10];
% XI = unique_heading(1):0.2:unique_heading(end);
% XI = [-10 -3.5:0.2:-1.3 -1.2:0.1:1.2 1.3:0.2:3.5 10];
XI = unique_heading';

if ~make_fake  % for now, only do all this if not making fake single-cue responses

for n = 1:cell_num 
    firing{1,2,1}(n,:) = interp1(X,resp_mean_all{n}{1,2,1},XI,method);
    firing{2,2,1}(n,:) = interp1(X,resp_mean_all{n}{2,2,1},XI,method);
    firing{2,2,2}(n,:) = interp1(X,resp_mean_all{n}{2,2,2},XI,method);
    firing{3,1,1}(n,:) = interp1(X,resp_mean_all{n}{3,1,1},XI,method);
    firing{3,2,1}(n,:) = interp1(X,resp_mean_all{n}{3,2,1},XI,method);
    firing{3,3,1}(n,:) = interp1(X,resp_mean_all{n}{3,3,1},XI,method);
    firing{3,1,2}(n,:) = interp1(X,resp_mean_all{n}{3,1,2},XI,method);
    firing{3,2,2}(n,:) = interp1(X,resp_mean_all{n}{3,2,2},XI,method);
    firing{3,3,2}(n,:) = interp1(X,resp_mean_all{n}{3,3,2},XI,method);
end

% % build covariance matrix r_noise to approximate a desired maximum noise
% % correlation, scaled by signal correlation
% for stim_type = 1:3   
%     if stim_type == 1
%         cohnum = 1;
%     else
%         cohnum = 2;
%     end
%     for i = 1:cell_num
%         for j = 1:cell_num
%             
%             %***** option 1: STIM-TYPE-SPECIFIC R_SIGNAL
% %             temp = corrcoef(resp_mean_all{i}{stim_type,2,cohnum}, resp_mean_all{j}{stim_type,2,cohnum});
%             %***** option 2: FIXED R_SIGNAL BASED ON (for example) VISUAL TUNING
% %             temp = corrcoef(resp_mean_all{i}{2,2,2}, resp_mean_all{j}{2,2,2});
%             %***** option 3: R_signal based on both vis and ves tuning, appended 
%             temp = corrcoef([resp_mean_all{i}{1,2,1} resp_mean_all{i}{2,2,2}], [resp_mean_all{j}{1,2,1} resp_mean_all{j}{2,2,2}]);            
%             r_signal{stim_type}(i,j) = temp(2,1);
%             if j==i
%                 r_noise{stim_type}(i,j) = 1; 
%             else
%                  r_noise{stim_type}(i,j) = noise_corr_max*exp( (r_signal{stim_type}(i,j)-1)/0.4785 );
%             end 
%         end
%     end
% end

for i = 1:cell_num
    for j = 1:cell_num
        %***** option 4: separate R_signal for vis and ves, but R_noise
        %(for all conditions) is based on mult. regress of both
        temp = corrcoef(resp_mean_all{i}{1,2,1}, resp_mean_all{j}{1,2,1});
        r_signal{1}(i,j) = temp(2,1);
        temp = corrcoef(resp_mean_all{i}{2,2,2}, resp_mean_all{j}{2,2,2});
        r_signal{2}(i,j) = temp(2,1);
        if j==i
            r_noise{1}(i,j) = 1;
            r_noise{2}(i,j) = 1;
            r_noise{3}(i,j) = 1;
        else
            r_noise{1}(i,j) = 0.01475 + 0.0855*r_signal{1}(i,j) + 0.07039*r_signal{2}(i,j); % these values are from Yong, probably trained animals -- DOUBLE CHECK IF USED AGAIN IN REVISION
            r_noise{2}(i,j) = r_noise{1}(i,j);
            r_noise{3}(i,j) = r_noise{1}(i,j);
            if noise_corr_max == 0
                r_noise{1}(i,j) = 0;
                r_noise{2}(i,j) = 0;
                r_noise{3}(i,j) = 0;
            end
        end
    end
end

% now generate partially correlated responses according to r_noise
if noisemodel == 'gausian'  % (this method only works for Gaussian noise)
    
    s=1; j=2; k=1; % vestib
    rootQ = sqrtm(r_noise{s});
    for i=1:length(XI)
        for r = 1:repetition
            z = normrnd(0,1,cell_num,1);
            y = rootQ * z;
            y = y.*sqrt(fano*firing{s,j,k}(:,i)) + firing{s,j,k}(:,i); % variance is fano*mean, and SD is sqrt of that
            y(find(y<0)) = 0; % no negative reponses
            firing_trial{s,j,k,i}(:,r) = y;
        end
    end

    s=2; j=2; % visual
    for k = 1:2 % coherence
        rootQ = sqrtm(r_noise{s});
        for i=1:length(XI)
            for r = 1:repetition
                z = normrnd(0,1,cell_num,1);
                y = rootQ * z;
                y = y.*sqrt(fano*firing{s,j,k}(:,i)) + firing{s,j,k}(:,i); % variance is fano*mean, and SD is sqrt of that
                y(find(y<0)) = 0; % no negative reponses
                firing_trial{s,j,k,i}(:,r) = y;
            end
        end
    end

    s=3; % combined
    for j = 1:3 % conflict angle
        for k = 1:2 % coherence
            rootQ = sqrtm(r_noise{s});
            for i=1:length(XI)
                for r = 1:repetition
                    z = normrnd(0,1,cell_num,1);
                    y = rootQ * z;
                    y = y.*sqrt(fano*firing{s,j,k}(:,i)) + firing{s,j,k}(:,i); % variance is fano*mean, and SD is sqrt of that
                    y(find(y<0)) = 0; % no negative reponses
                    firing_trial{s,j,k,i}(:,r) = y;
                end
            end
        end
    end

else % for Poisson, just generate firing_trial (zero correlation) using poissrnd at the interpolated values
    
    for n = 1:cell_num 
        for s = 1:3; % stim type
            for j = 1:3 % conflict angle
                for k = 1:2 % coherence
                    if ~isempty(firing{s,j,k})
                        for i=1:length(XI)
                            firing_trial{s,j,k,i}(n,1:repetition) = poissrnd(firing{s,j,k}(n,i),1,repetition);
                        end
                    end
                end
            end
        end
    end
    
end

end

figure; contourf(r_noise{3}); caxis([-0.2 0.2]); colorbar;
cormat = r_noise{3};