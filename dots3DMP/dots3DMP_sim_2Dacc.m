% simulate "CCC" (cue-combination + confidence) experiment
% specifically in the case of the vis-ves heading task

% CF started it 2016
% heavily modified early 2019
% switched to 2D accumulator model 06/2020

% build expt and hand-pick some model params

clear
close all

RTtask = 1;
conftask = 2; % 1 - sacc endpoint, 2 - PDW

plotExampleTrials = 1;

nreps = 200; % number of repetitions of each unique trial type
            % start small to verify it's working, then increase
            % (ntrials depends on num unique trial types)

cohs = [0.3 0.7]; % visual coherence levels (these are really just labels, since k's are set manually)
hdgs = [-10 -5 -2.5 -1.25 0 1.25 2.5 5 10]; % heading angles
% hdgs = [-10 -3.5 -1.25 1.25 3.5 10]; % heading angles
deltas = [-3 0 3]; % conflict angle; positive means vis to the right
mods = [1 2 3]; % stimulus modalities: ves, vis, comb
% duration = 2000; % stimulus duration (ms)
duration = 1300; 

% sensitivity constants for converting heading angle into mean of momentary evidence
ks = 25; % scale factor, for quickly testing different k levels
  % set manually to get reasonable results from images_dtb_2d
kves = ks;
kvis = ks * [2 4.5]/3; % straddles vestibular reliability, by construction

theta = 0.7; % threshold for high bet in logOdds, ignored if conftask==1
    
sigmaVes = 0.01; % std of momentary evidence
sigmaVis = [0.01 0.01]; % allow for separate sigmas for condition, coherence
    % set manually to get reasonable looking dv trajectories;
    % also affects peakiness/flatness of confidence curve
    
B = 1.2; % bound height (also hand-tuned)

% draw non-decision times from Gaussian dist
% muTnd will be a parameter to fit, sdTnd can be fixed at reasonable value
muTnd = 300; sdTnd = 60; % ms
% Tnd = 0.4; % fixed val

% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis]);

R.t = 0.001:0.001:duration/1000;
R.Bup = B;
R.drift = k * sind(hdgs(hdgs>=0)); % takes only unsigned drift rates
R.lose_flag = 1;
R.plotflag = 0; % 1 = plot, 2 = plot and export_fig

P =  images_dtb_2d(R);
% uses method of images to compute PDF of the 2D DV, as in van den Berg et
% al. 2016 (similar to Kiani et al. 2014)

% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
% 07/2020 lab settings...160cm in 1.3s, sigma=0.14

ampl = 0.16; % movement in metres
pos = normcdf(1:duration,duration/2,0.14*duration)*ampl;
vel = gradient(pos)*1000; % metres/s
acc = gradient(vel);

% vel = normpdf(1:duration,duration/2,210);
% vel = 0.37*vel./max(vel);
% acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)

origParams.ks = ks;
origParams.sigma = sigmaVes;
origParams.B = B;
origParams.Tnd = muTnd;
if conftask==2
    origParams.theta = theta; % PDW only
end

%% build trial list

[hdg, modality, coh, delta, ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps);

% sample durations from truncated exponential?
    % not practical experimentally.
    % instead, make constant dur, assume RT task, or (equivalently, as far as
    % the simulation is concerned) fixed duration with internal bounded
    % process (Kiani et al. 2008)
dur = ones(ntrials,1) * duration;

Tnds = muTnd + randn(ntrials,1).*sdTnd;

%% bounded evidence accumulation

% assume momentary evidence is proportional to sin(heading),
% as in drugowitsch et al 2014

% dv_all = cell(ntrials,1); % shouldn't need to store every trial's DV, but if you want to, it's here

choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1); % now this is the value of the losing accumulator
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

% for plotting
doneWith1 = 0;
doneWith2 = 0;
doneWith3 = 0;

% ME is now a draw from bivariate normal with mean vector Mu and covariance matrix V
% start with correlation matrix:
S = [1 -1/sqrt(2) ; -1/sqrt(2) 1];
% -1/sqrt(2) is the correlation for our version of images_dtb;
% can only be changed with an update to the flux file

tic
for n = 1:ntrials
    Tnd = Tnds(n) / 1000; % Tnd for nth trial in seconds
    
    switch modality(n)
        case 1
            mu = acc .* kves * sind(hdg(n)) / 1000; % mean of momentary evidence
                % (I'm guessing drift rate in images_dtb is per second, hence div by 1000)
            s = [sigmaVes sigmaVes]; % standard deviaton vector (see below)
        case 2
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n)) / 1000;
            s = [sigmaVis(cohs==coh(n)) sigmaVis(cohs==coh(n))];
        case 3
            % positive delta defined as ves to the left, vis to the right
            muVes = acc .* kves               * sind(hdg(n)-delta(n)/2) / 1000;
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n)+delta(n)/2) / 1000;

            % optimal weights (Drugo et al.) 
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            mu = wVes.*muVes + wVis.*muVis;

            % the DV is a sample from a dist with mean = weighted sum of
            % means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
            s = [sigmaComb sigmaComb];
    end

%     Mu = [mu,-mu]; % mean vector for 2D DV
    Mu = [mu; -mu]';
    
    % convert correlation to covariance matrix
    V = diag(s)*S*diag(s);
%     dv = [0 0; cumsum(mvnrnd(Mu,V,dur(n)-1))]; % bivariate normrnd
    dv = [0 0; cumsum(mvnrnd(Mu,V))]; % dv is now scaled by physical signal vel/acc
    
    % because Mu is signed according to heading (positive=right),
    % dv(:,1) corresponds to evidence favoring rightward, not evidence
    % favoring the correct decision (as in Kiani eqn. 3 and images_dtb)

%     dv_all{n} = dv;
    % decision outcome
    cRT1 = find(dv(:,1)>=B, 1);
    cRT2 = find(dv(:,2)>=B, 1);
    % the options are:
    % (1) only right accumulator hits bound,
    if ~isempty(cRT1) && isempty(cRT2)
        RT(n) = cRT1/1000 + Tnd;
        finalV(n) = dv(cRT1,2); % only 1 hit, so 2 is the loser
        hitBound(n) = 1;
        choice(n) = 1;
    % (2) only left accumulator hits bound, 
    elseif isempty(cRT1) && ~isempty(cRT2)
        RT(n) = cRT2/1000 + Tnd;
        finalV(n) = dv(cRT2,1); % only 2 hit, so 1 is the loser
        hitBound(n) = 1;
        choice(n) = -1;
    % (3) neither hits bound,
    elseif isempty(cRT1) && isempty(cRT2)
        RT(n) = dur(n)/1000 + Tnd;
            % this is interesting: which DV matters for confidence
            % if neither hits bound? for now take the (abs) maximum,
            % ie whichever was closest to hitting bound. (alternative
            % would be their average?)
            % SJ 07/2020 finalV is fully determined by distance of loser from bound when winner
            % hits, so in this case, should also account for where winner is wrt bound 
            % imagine winner 'did' hit bound, then where
            % would loser be relatively speaking (since logOdds map is
            % fixed)
           
        whichWon = dv(end,:)==max(dv(end,:));
        finalV(n) = dv(end,~whichWon) + B-dv(end,whichWon);
        % ^ effectively shifting the losing dv up by whatever the
        % difference is between the bound and the winning dv
%         finalV(n) = dv(end,~whichWon); % the not-whichWon is the loser
        % % finalV(n) = mean(dvEnds); 
        
        hitBound(n) = 0;
        a = [1 -1];
        choice(n) = a(whichWon);
    % (4) or both do
    else
        RT(n) = min([cRT1 cRT2])/1000 + Tnd;
        whichWon = [cRT1<=cRT2 cRT1>cRT2];
        finalV(n) = dv(min([cRT1 cRT2]),~whichWon); % the not-whichWon is the loser
        hitBound(n) = 1;
        a = [1 -1];
        choice(n) = a(whichWon);
    end
        
    % use map to look up log-odds that the motion is rightward
    diffV = abs((P.y+B)-finalV(n));
    diffT = abs(R.t-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = P.logOddsCorrMap(thisV(1), thisT(1));
    
    if conftask==1 % sacc endpoint
        expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
        conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1
    elseif conftask==2 % PDW
        conf(n) = logOddsCorr(n) > theta;
    end
    if isnan(conf(n)), conf(n)=0; end % if dvs are almost overlapping, force conf to zero as it can sometimes come out as NaN
    
    if plotExampleTrials

        if modality(n)==1 && hdg(n)==5 && choice(n)==1 && doneWith1==0 % make a better plot for talk/poster;
                              % must not shuffle trial list for this to go in correct in order
            
            figure(1000); set(gcf, 'Color', [1 1 1], 'Position', [100 100 350 375], 'PaperPositionMode', 'auto'); clf;
            hold on; box off;
            xlabel('Time (ms)');
            ylabel('Accum. evidence');
            ylim([-1.25 1.1].*B); xlim([0 duration]);
            %             set(gca,'yTick',-0.5:0.5:2);
            set(gca,'xTick',0:500:2000);
            changeAxesFontSize(gca,18,18);
            
            yl = get(gca,'ylim');
            plot(1:length(vel),vel+yl(1),'k--','linew',0.5)
            
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'k-','LineWidth',2); hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'color',[1 1 1]*0.5,'LineWidth',1); hold on; % losing (L)
            fprintf('Vest hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));
        
            plot(1:length(dv),ones(1,length(dv))*B,'k-','LineWidth',4);
    
            plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'k:')
%             plot([round(RT(n)*1000) round(RT(n)*1000)],[dv(round((RT(n)-Tnd)*1000),2) B],'r--')

%             text((RT(n)-Tnd/2)*1000,1,'NDT','horizo','center','verti','bottom','fontweight','bold','fontsize',14)
%             keyboard
%             text((RT(n)+0.1-Tnd)*1000,0,sprintf('Conf = %.2f',conf(n)),'fontweight','bold','fontsize',14)
%             text(RT(n)*1000-Tnd,1.12*B,sprintf('RT = %.2fs',RT(n)),'fontweight','bold','fontsize',14)
            
            box off;
            xlabel('Time (ms)');
%             ylabel('Accum. evidence for leftward');
            %ylim([-0.4*B B*1.1]); xlim([0 duration]);
%             set(gca,'yTick',-0.5:0.5:2);
            %set(gca,'xTick',0:500:2000);
            %changeAxesFontSize(gca,18,18);
            %export_fig('ves_acc2','-eps');

            doneWith1=1; n
        end
        if modality(n)==2 && hdg(n)==5 && coh(n)==cohs(2) && choice(n)==1 && doneWith2==0
            figure(1000);
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'r-','LineWidth',2); %hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'m-','LineWidth',1); % losing (L)            doneWith2 = 1;
            plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'r:')
            fprintf('Vis hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));
            doneWith2 = 1; n
        end
        if modality(n)==3 && hdg(n)==5 && coh(n)==cohs(2) && choice(n)==1 && delta(n)==0 && doneWith3==0
            figure(1000);
            plot(dv(1:round((RT(n)-Tnd)*1000),1),'b-','LineWidth',2); %hold on; % winning (R, since we preselect for R choice)
            plot(dv(1:round((RT(n)-Tnd)*1000),2),'c-','LineWidth',1); % losing (L)
            plot([round((RT(n)-Tnd)*1000) round((RT(n)-Tnd)*1000)],dv(round((RT(n)-Tnd)*1000),:),'b:')
            fprintf('Comb hdg = %.1f, RT (-NDT) = %.2f, conf = %.2f',hdg(n),RT(n)-Tnd,conf(n));
            doneWith3 = 1; n
%             export_fig('threeTrials','-eps');
        end
    end

end
toc

figure(1000);
plot(1:length(vel),vel+yl(1),'k--','linew',0.5)
choice(choice==0) = sign(randn); % not needed under usual circumstances

% sanity check:
pCorrect_total = (sum(choice==1 & hdg>0) + sum(choice==-1 & hdg<0)) / ntrials





%% format data like the real expt

choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

data.modality = modality;
data.heading = hdg;
data.coherence = coh;
data.delta = delta;
data.choice = choice;
data.RT = RT; % already in seconds
data.conf = conf;

if conftask==2
    data.PDW=data.conf;
end

% cd('/Users/stevenjerjian/Desktop/FetschLab/DATA analysis')
% save('2DAcc_simdata_sacc_stimprof.mat','data','cohs','deltas','hdgs','mods','origParams','RTtask','conftask')

%% plots

dots3DMP_parseData
dots3DMP_plots

% dots3DMP_parseData_splitConf
% dots3DMP_plots_splitConf

if 0
%% fit cumulative gaussians

dots3DMP_fit_cgauss
dots3DMP_plots_cgauss
% dots3DMP_plots_cgauss_forTalk % nicer looking versions


%% now try fitting the fake data to recover the generative parameters

options.fitMethod = 'fms';
% options.fitMethod = 'global';
% options.fitMethod = 'multi';
% options.fitMethod = 'pattern';
% options.fitMethod = 'bads';

% initial guess (or hand-tuned params)
ks      = 22;
sigma   = 0.01;
B       = 1.2;
Tnd     = 300;
fixed   = [0 1 1 1];
guess   = [ks sigma B Tnd];

if conftask==2 % PDW
    theta = 1; 
    
    fixed = [0 1 0 1 1];
    guess = [ks sigma B Tnd theta];
end

% ************************************
% set all fixed to 1 for hand-tuning:
% fixed(:)=1;
% (can be used to fix some params and not others)
% ************************************

% plot error trajectory (prob doesn't work with parallel fit methods)
options.ploterr = 0;

options.RTtask = RTtask;
options.conftask = conftask; % 1 - sacc endpoint, 2 - PDW

[X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed);

% plot it!
% dots3DMP_plots_fit(data,fitInterp,conftask,RTtask)

end



