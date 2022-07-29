% Standardized computation of marginals and logOddsCorrect map for 1D DDM
% (FP4.mex, Kiani et al. 2009 method)

% Adapted from errfcn_DDM_1D_wConf.m, replaces makeLogOddsCorrMap_smooth
% for use in simulations, so that fitting and sims use identical code, e.g.
% for parameter recovery.

% Script for now, could convert to a function

% CF 07/2022


if ~exist('coh','var'); coh = data.scoh; end
if ~exist('max_dur','var'); max_dur = max(data.dur); end


% get the stimulus coh levels, the frequency of each level and the maximum stimulus duration
% Note that coh must be signed for the algorithm to work properly 
coh_set = unique(coh);
coh_set_freq = nan(length(coh_set),1);        
for c = 1 : length(coh_set)
    coh_set_freq(c) = sum(coh==coh_set(c))/length(coh);
end

    % define the grid resolution for time and decision variable
dt = 0.1; % 0.1 ms seems to work best for FP4, even though we don't store the vars at this resolution
dx = min(0.1, B/100);
    % define the time axis 
t = dt:dt:max_dur;
    % define xmesh and the ghost cells for bound crossing probability
b_margin = repmat(4*sigma, [1 2]);
xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';
vAxis = xmesh; % legacy

    % adjust dx so that the mesh has zero and is therefore symmetric 
while ~any(xmesh==0)
    [~, I] = min(abs(xmesh));
    delta_mesh = xmesh(I);
    b_margin = b_margin + delta_mesh;
    xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)'; % recompute the mesh
end
if mod(length(xmesh),2)~=1
    error('the length of xmesh must be odd');
end

    % define a delta function on xmesh
delta = zeros(size(xmesh));
delta(abs(xmesh)==min(abs(xmesh))) = 1;

    % initialize:
% probability of bound crossing as a function of time
Ptb_coh = zeros(max_dur, 2, length(coh_set));                % time * bound * signed_coh
% % Ptb_coh_forMarginals = Ptb_coh;
% probability density of decision variable (x) across time
Pxt_coh = zeros(length(xmesh), max_dur, length(coh_set));    % xmesh * time * signed_coh
% % Pxt_coh_forMarginals = Pxt_coh;
% marginal densities
Ptb_marginal = zeros(max_dur, 2, 2);                            % time * bound * motion_direction(1 left, 2 right)
Pxt_marginal = zeros(length(xmesh), max_dur, 2);                % xmesh * time * motion_direction
    % NOTE: time in these vars is at 1 ms resolution, not 0.1 (downsampled, see below)
    % NOTE: CF changed motion direction indices to match lower/upper bound indices (fixed output of FP4)

    %run FP4 to calculate Ptb and Pxt for each coherence and stimulation condition 
for c = 1 : length(coh_set)
    % All arguments to FP4 are corrupted by FP4 (allegedly). Reset them w/in the loop, before each call
    X = xmesh;
    uinit = delta; % start with delta function
    Mu = k * coh_set(c);
    Sig = sigma;
    [b_change, ~] = bcinit(B,t,dt,'flat',NaN,NaN); % b_change is for time-varying (e.g. collapsing) bounds
    b_marg = b_margin;
    dT = dt;
    
    %%%%
    [~, ~, Ptb, ~, Pxt] = FP4(X, uinit, Mu, Sig, b_change, b_marg, dT);
    %%%%
    
        %store the arrays for each coherence level, use 1 ms time resolution (downsample)
    Ptb_coh(:,:,c) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
    Pxt_coh(:,:,c) = Pxt(:,1/dt:1/dt:end); % this is downsampled but Ptb is local-summed! what gives??!?!

%         %store the arrays for each coherence level (if dt=1ms)
%     Ptb_coh(:,:,c) = [Ptb(2:end,1) Ptb(2:end,1)];
%     Pxt_coh(:,:,c) = Pxt;

    % store separate vars for the marginals, to allow them to be computed
    % differently ie for more complicated models, or to downsample further;
% %     uinit = delta;      %start with delta function
% %     [b_change, ~] = bcinit(B,t,dt,'flat',NaN,NaN);
% %     [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin, dt);
    % currently we just use the same distributions, so save memory and
    % don't duplicate them here
% % %     Ptb_coh_forMarginals(:,:,c) = Ptb_coh(:,:,c);
% % %     Pxt_coh_forMarginals(:,:,c) = Pxt_coh(:,:,c); 
end

    % calculate the marginals (marginalize over coherence, separately for
    % each motion direction -- this depends on frequency of presentation)
for c = 1 : length(coh_set)
    F = coh_set_freq(coh_set==coh_set(c));
        %now use Ptb and Pxt to calculate/update the marginal. 
%     Pxt = Pxt_coh_forMarginals(:,:,c);
%     Ptb = Ptb_coh_forMarginals(:,:,c);
    Pxt = Pxt_coh(:,:,c); % see above; no need to duplicate
    Ptb = Ptb_coh(:,:,c);
        
    if coh_set(c)>0         %rightward motion (index 2 is upper bound = right)
        Ptb_marginal(:,:,2) = Ptb_marginal(:,:,2) + F*Ptb; % time * bound
        Pxt_marginal(:,:,2) = Pxt_marginal(:,:,2) + F*Pxt; % xmesh * time
    elseif coh_set(c)<0     %leftward motion (index 1 is lower bound = left)
        Ptb_marginal(:,:,1) = Ptb_marginal(:,:,1) + F*Ptb; % time * bound
        Pxt_marginal(:,:,1) = Pxt_marginal(:,:,1) + F*Pxt; % xmesh * time
    else                    %ambiguous motion
        Ptb_marginal(:,:,2) = Ptb_marginal(:,:,2) + 0.5*F*Ptb;
        Pxt_marginal(:,:,2) = Pxt_marginal(:,:,2) + 0.5*F*Pxt;
        Ptb_marginal(:,:,1) = Ptb_marginal(:,:,1) + 0.5*F*Ptb;
        Pxt_marginal(:,:,1) = Pxt_marginal(:,:,1) + 0.5*F*Pxt;
    end
end

    % find out which combination of DV and time is associated with high/low wager based on theta
%if bound crossing does not happen
bet_high_xt = nan(length(xmesh), max_dur);

I = xmesh>0;
logPosteriorOddsRight = log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1));
bet_high_xt(I,:) = logPosteriorOddsRight > theta;

I = xmesh<0;
logPosteriorOddsLeft = log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2));
bet_high_xt(I,:) = logPosteriorOddsLeft > theta2;

I = xmesh==0; % dv=0, assume a 50/50 guess when that happens
if sum(I)>0
    bet_high_xt(I,:) = 0>(theta+theta2)/2;
    LPOmiddle = mean([log(Pxt_marginal(find(I)+1,:,1)./Pxt_marginal(find(I)+1,:,2)) ; log(Pxt_marginal(find(I)-1,:,1)./Pxt_marginal(find(I)-1,:,2))]);
    logOddsCorrMap = [logPosteriorOddsLeft ; LPOmiddle ; logPosteriorOddsRight];
else
    logOddsCorrMap = [logPosteriorOddsLeft ; logPosteriorOddsRight];
end

%if bound crossing happens 
bet_high_tb = nan(max_dur, 2);
                     % Ptb_marginal: time * bound * motion_direction (1 left, 2 right)
bet_high_tb(:,1) = log(Ptb_marginal(:,1,1)./Ptb_marginal(:,1,2)) > theta2;    %lower bound crossing (left)
bet_high_tb(:,2) = log(Ptb_marginal(:,2,2)./Ptb_marginal(:,2,1)) > theta;     %upper bound crossing (right)

TBint1 = find(bet_high_tb(:,1)==0,1,'first'); % theta-bound intersection point
TBint2 = find(bet_high_tb(:,2)==0,1,'first');
if ~isempty(TBint1); TBint(1) = TBint1; else TBint(1) = NaN; end
if ~isempty(TBint2); TBint(2) = TBint2; else TBint(2) = NaN; end
    % WARNING: these intersections don't quite match up to where
    % logOddsCorrMap > theta just below the bound (seems not to matter??)

% plot the marginal PDFs, logOddsCorr, and high/low bet regions
if options.plot
    n=30; % number of color levels (smoothness)
    figure(123); set(gcf,'Color',[1 1 1],'Position',[731 568 934 766],'PaperPositionMode','auto'); clf;
    x = repmat(1:size(logOddsCorrMap,2),size(logOddsCorrMap,1),1)';
    y = repmat(xmesh,1,size(logOddsCorrMap,2))';

    subplot(2,2,1);         
    temp = Pxt_marginal(:,:,2); temp(temp<1e-10) = 1e-10;
    z = log10(temp)';
    [~,h] = contourf(x,y,z,n); title('log(PxtMarg), rightward'); colorbar;
    if n>20; set(h,'LineColor','none'); end
    xlabel('t (ms)'); ylabel('x');
    xlim([1 1000]);

    subplot(2,2,2);
    temp = Pxt_marginal(:,:,1); temp(temp<1e-10) = 1e-10;
    z = log10(temp)';
    [~,h] = contourf(x,y,z,n); title('log(PxtMarg), leftward'); colorbar;
    if n>20; set(h,'LineColor','none'); end
    xlabel('t (ms)'); ylabel('x');
    xlim([1 1000]);

    subplot(2,2,3);
    z = logOddsCorrMap';
    [~,h] = contourf(x,y,z,n); colorbar; title('Log odds correct');
    if n>20; set(h,'LineColor','none'); end
    xlabel('t (ms)'); ylabel('x');
    xlim([1 1000]);

    subplot(2,2,4);
    z = bet_high_xt';
    [~,h] = contourf(x,y,z,n); colorbar; title('bet high');
    if n>20; set(h,'LineColor','none'); end
    xlabel('t (ms)'); ylabel('x');
    xlim([1 1000]);
end
