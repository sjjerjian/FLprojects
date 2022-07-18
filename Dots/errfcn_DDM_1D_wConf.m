function [err,fit,parsedFit] = errfcn_DDM_1D_wConf(param, guess, fixed, data, options)
    
tic

global call_num

param = getParam(param, guess, fixed);

k = param(1);
B = abs(param(2)); % don't accept negative bound heights
sigma = 1;
theta = abs(param(3)); %or negative thetas
theta2 = theta; % allows expanded model with separate thetas for each choice (currently unused)
alpha = param(4); % base rate of low-conf choices

% get the stimulus coh levels, the frequency of each level and the maximum stimulus duration
% Note that coh must be signed for the algorithm to work properly 
coh_set = unique(data.scoh);
coh_set_freq = nan(length(coh_set),1);        
for s = 1 : length(coh_set)
    coh_set_freq(s) = sum(data.scoh==coh_set(s))/length(data.scoh);
end

%     % Miguel's idea: freq for marg should be based on unsigned coh
%     coh_set_unsigned = unique(abs(data.scoh));
%     coh_set_freq = nan(length(coh_set_unsigned),1);        
%     for s = 1 : length(coh_set_unsigned)
%         coh_set_freq(s) = sum(abs(data.scoh)==coh_set_unsigned(s))/length(data.scoh);
%     end

    % define the grid resolution for time and decision variable
dt = 0.1; % 0.1 ms seems to work best for FP4, even though we don't store the vars at this resolution
dx = min(0.1, B/100);
    % define the time axis 
max_dur = max(data.dur);
t = dt:dt:max_dur;
    % define xmesh and the ghost cells for bound crossing probability
b_margin = repmat(4*sigma, [1 2]);
xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';

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
Ptb_coh_forMarginals = Ptb_coh;
% probability density of decision variable (x) across time
Pxt_coh = zeros(length(xmesh), max_dur, length(coh_set));    % xmesh * time * signed_coh
Pxt_coh_forMarginals = Pxt_coh;
% marginal densities
Ptb_marginal = zeros(max_dur, 2, 2);                            % time * bound * motion_direction(1 left, 2 right)
Pxt_marginal = zeros(length(xmesh), max_dur, 2);                % xmesh * time * motion_direction
    % NOTE t axis is at 1 ms resolution, not dt (downsampled, see below)

    % NOTE: CF changed motion direction indices to match lower/upper
    % bound indices (fixed output of FP4)

    %run FP4 to calculate Ptb and Pxt for each coherence and stimulation condition 
for s = 1 : length(coh_set)
    mu = k * coh_set(s);
    uinit = delta; % start with delta function
    % b_change is for time-varying (e.g. collapsing) bounds, which must be re-initialized after each call of FP4
    [b_change, ~] = bcinit(B,t,dt,'flat',NaN,NaN);
    [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin, dt);
        %store the arrays for each coherence level, use 1 ms time resolution (downsample)
    Ptb_coh(:,:,s) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
    Pxt_coh(:,:,s) = Pxt(:,1/dt:1/dt:end);
end

    % repeat for marginals [in a separate loop to allow them to be computed
    % differently, for some more complicated models]
for s = 1 : length(coh_set)
    mu = k * coh_set(s);
    uinit = delta;      %start with delta function
    [b_change, ~] = bcinit(B,t,dt,'flat',NaN,NaN);
    [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin, dt);
    Ptb_coh_forMarginals(:,:,s) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
    Pxt_coh_forMarginals(:,:,s) = Pxt(:,1/dt:1/dt:end);
end

    % calculate the marginals (marginalize over coherence, separately for
    % each motion direction -- this depends on frequency of presentation)
for s = 1 : length(coh_set)
    F = coh_set_freq(coh_set==coh_set(s));
% % %         F = coh_set_freq(coh_set_unsigned==abs(coh_set(s))); % MVL
        %now use Ptb and Pxt to calculate/update the marginal. 
    Pxt = Pxt_coh_forMarginals(:,:,s);
    Ptb = Ptb_coh_forMarginals(:,:,s);
    if coh_set(s)>0         %rightward motion (index 2 is upper bound = right)
        Ptb_marginal(:,:,2) = Ptb_marginal(:,:,2) + F*Ptb; % time * bound
        Pxt_marginal(:,:,2) = Pxt_marginal(:,:,2) + F*Pxt; % xmesh * time
    elseif coh_set(s)<0     %leftward motion (index 1 is lower bound = left)
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
    plotMiddle = mean([log(Pxt_marginal(find(I)+1,:,1)./Pxt_marginal(find(I)+1,:,2)) ; log(Pxt_marginal(find(I)-1,:,1)./Pxt_marginal(find(I)-1,:,2))]);
    plotAll = [logPosteriorOddsLeft ; plotMiddle ; logPosteriorOddsRight];
else
    plotAll = [logPosteriorOddsLeft ; logPosteriorOddsRight];
end

%if bound crossing happens 
bet_high_tb = nan(max_dur, 2);
                       % Ptb_marginal: time * bound * motion_direction (1 right, 2 left)
%     bet_high_tb(:,1) = log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1)) > theta2;    %lower bound crossing
%     bet_high_tb(:,2) = log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2)) > theta;     %upper bound crossing
bet_high_tb(:,1) = log(Ptb_marginal(:,1,1)./Ptb_marginal(:,1,2)) > theta2;    %lower bound crossing
bet_high_tb(:,2) = log(Ptb_marginal(:,2,2)./Ptb_marginal(:,2,1)) > theta;     %upper bound crossing

thetaIntersectionPoint(1) = find(bet_high_tb(:,1)==0,1,'first');
thetaIntersectionPoint(2) = find(bet_high_tb(:,2)==0,1,'first');


% plot the marginals, logOddsCorr, and high/low bet regions
if options.plot
    n=30; % number of color levels (smoothness)
    figure(123); set(gcf,'Color',[1 1 1],'Position',[731 568 934 766],'PaperPositionMode','auto'); clf;
    x = repmat(1:size(plotAll,2),size(plotAll,1),1)';
    y = repmat(xmesh,1,size(plotAll,2))';

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
    z = plotAll';
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


% accumulate log likelihood
expectedPright = nan(size(data.scoh));
expectedPrightHigh = nan(size(data.scoh));
expectedPrightLow = nan(size(data.scoh));
expectedPhigh = nan(size(data.scoh));
expectedPhighCorr = nan(size(data.scoh));
expectedPhighErr = nan(size(data.scoh));
expectedRT = nan(size(data.scoh));
expectedRThigh = nan(size(data.scoh));
expectedRTlow = nan(size(data.scoh));
t_ms = 1:size(Ptb,1);
pRight_model = nan(length(coh_set),1);
pHigh_model = nan(length(coh_set),1);
pHighCorr_model = nan(length(coh_set),1);
pHighErr_model = nan(length(coh_set),1);
pRightHigh_model = nan(length(coh_set),1);
pRightLow_model = nan(length(coh_set),1);
meanRT_model = nan(length(coh_set),1);
meanRThigh_model = nan(length(coh_set),1);
meanRTlow_model = nan(length(coh_set),1);
try
    n = 0;
    LL = 0;
    for s = 1 : length(coh_set)

        Ptb = Ptb_coh(:,:,s); % time * bound
        Pxt = Pxt_coh(:,:,s); % xmesh * time
        
        % calculate probabilities of the four outcomes: right/left x high/low,
        % as a function of time (dur)
            % *** I don't think this makes sense for RT (periDW): Pxt
            % shouldn't enter into it, because all trials end at the bound
        PrightHigh = cumsum(Ptb(:,2).*(bet_high_tb(:,2)==1)) + ...
                    sum(Pxt(xmesh>0,:).*(bet_high_xt(xmesh>0,:)==1),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==1))'; % half the probability of dv=0 (assuming a 50/50 guess when that happens)
        PrightLow = cumsum(Ptb(:,2).*(bet_high_tb(:,2)==0)) + ...
                    sum(Pxt(xmesh>0,:).*(bet_high_xt(xmesh>0,:)==0),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==0))';
        PleftHigh = cumsum(Ptb(:,1).*(bet_high_tb(:,1)==1)) + ...
                    sum(Pxt(xmesh<0,:).*(bet_high_xt(xmesh<0,:)==1),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==1))';
        PleftLow =  cumsum(Ptb(:,1).*(bet_high_tb(:,1)==0)) + ...
                    sum(Pxt(xmesh<0,:).*(bet_high_xt(xmesh<0,:)==0),1)' + ...
                    0.5*(Pxt(xmesh==0,:).*(bet_high_xt(xmesh==0,:)==0))';

% % %             % temp: plot the probability of the 4 outcomes as a function of time
% % %             figure; hold on; 
% % %             plot(PrightHigh,'b-'); plot(PrightLow,'b--');
% % %             plot(PleftHigh,'r-'); plot(PleftLow,'r--');
% % %             legend('R high','R low','L high','L low');
% % %             figure; hold on; 
% % %             plot(PrightHigh+PleftHigh,'g-');
% % %             plot(PrightLow+PleftLow,'m-');
% % %             legend('HIGH','LOW');
% % %             figure; hold on; 
% % %             plot(PrightHigh+PrightLow,'c-');
% % %             plot(PleftHigh+PleftLow,'k-');
% % %             legend('RIGHT','LEFT');


        %ensure total prob sums to one
        Ptot = PrightHigh + PrightLow + PleftHigh + PleftLow;
        if any(abs(Ptot-1)>1e-3)
            warning('the probabilities (Pright, Pleft, and Pstrg) do not add up to one!\n\tPtot=%f (for s:%d)\n', min(Ptot), s);
            PrightHigh = PrightHigh./Ptot;
            PrightLow = PrightLow./Ptot;
            PleftHigh = PleftHigh./Ptot;
            PleftLow = PleftLow./Ptot;
        end

        % adjust the probabilities for the base rate of low-conf bets:
        % the idea is that Phigh and Plow each get adjusted down/up in
        % proportion to how close they are to 1 or 0, respectively
        Phigh = PrightHigh + PleftHigh;
        Plow = PrightLow + PleftLow;
        Phigh2 = Phigh - alpha*Phigh;
        Plow2 = Plow + alpha*(1-Plow);

        % now recover the right/left proportions
        PrightHigh = (PrightHigh./Phigh) .* Phigh2;
        PrightLow = (PrightLow./Plow) .* Plow2;
        PleftHigh = (PleftHigh./Phigh) .* Phigh2;
        PleftLow = (PleftLow./Plow) .* Plow2;

        % offset to avoid negative P, then renormalize
        % [this should no longer be necessary!]
        minP = min([PrightHigh ; PrightLow ; PleftHigh; PleftLow]);
        if minP<0
            warning('shouldn''t happen'); keyboard
            PrightHigh = PrightHigh - minP;
            PrightLow = PrightLow - minP; 
            PleftHigh = PleftHigh - minP;
            PleftLow = PleftLow - minP;
        end
        Ptot = PrightHigh + PrightLow + PleftHigh + PleftLow;
        PrightHigh = PrightHigh./Ptot;
        PrightLow = PrightLow./Ptot;
        PleftHigh = PleftHigh./Ptot;
        PleftLow = PleftLow./Ptot;

        % calculate expected p(Right) and P(high), given the durs
            % *** not sure this can be right for RT version, see above
        I = data.scoh==coh_set(s);
        dur = data.dur(I);
        expectedPright(I) = (PrightHigh(dur) + PrightLow(dur)) ./ (PrightHigh(dur) + PrightLow(dur) + PleftHigh(dur) + PleftLow(dur));
            pRight_model(s) = mean(expectedPright(I)); % this feels wrong, somehow
        expectedPhigh(I) = (PrightHigh(dur) + PleftHigh(dur)) ./ (PrightHigh(dur) + PrightLow(dur) + PleftHigh(dur) + PleftLow(dur));
            pHigh_model(s) = mean(expectedPhigh(I));
        expectedPrightHigh(I) = PrightHigh(dur) ./ (PrightHigh(dur) + PleftHigh(dur)); % NOT simply PrightHigh as a quaternary outcome, because it's Pright *conditioned* on high
            pRightHigh_model(s) = mean(expectedPrightHigh(I));
        expectedPrightLow(I) = PrightLow(dur) ./ (PrightLow(dur) + PleftLow(dur));
            pRightLow_model(s) = mean(expectedPrightLow(I));
                        
        % new: pHigh separated by correct/error
        if coh_set(s)>=0 % rightward trials
            expectedPhighCorr(I) = PrightHigh(dur) ./ (PrightHigh(dur) + PrightLow(dur));
            expectedPhighErr(I) = PleftHigh(dur) ./ (PleftHigh(dur) + PleftLow(dur));
        else % leftward trials
            expectedPhighCorr(I) = PleftHigh(dur) ./ (PleftHigh(dur) + PleftLow(dur));
            expectedPhighErr(I) = PrightHigh(dur) ./ (PrightHigh(dur) + PrightLow(dur));
        end            
        pHighCorr_model(s) = mean(expectedPhighCorr(I));
        pHighErr_model(s) = mean(expectedPhighErr(I));
            
        % compute mean RT, just for consistency w 2D (not used for fitting;
        % effectively this is a fit-then-predict strategy for RT tasks)
        cohI = double(coh_set(s)>=0)+1;
        expectedRT(I) = (sum(Ptb(:,cohI) .* t_ms') / sum(Ptb(:,cohI)))/1000 + 0.35; % change to s and add a reasonable Tnd
            meanRT_model(s) = mean(expectedRT(I));
        expectedRThigh(I) = (sum(Ptb(1:thetaIntersectionPoint(cohI)-1,cohI) .* t_ms(1:thetaIntersectionPoint(cohI)-1)') / sum(Ptb(1:thetaIntersectionPoint(cohI)-1,cohI)))/1000 + 0.35;
            meanRThigh_model(s) = mean(expectedRThigh(I));
        expectedRTlow(I) = (sum(Ptb(thetaIntersectionPoint(cohI):end,cohI) .* t_ms(thetaIntersectionPoint(cohI):end)') / sum(Ptb(thetaIntersectionPoint(cohI):end,cohI)))/1000 + 0.35;
            meanRTlow_model(s) = mean(expectedRTlow(I));
                
        % update log-likelihood
        dur_rightHigh = data.dur(I & data.choice==1 & data.PDW==1);
        dur_rightLow = data.dur(I & data.choice==1 & data.PDW==0);
        dur_leftHigh = data.dur(I & data.choice==0 & data.PDW==1);
        dur_leftLow = data.dur(I & data.choice==0 & data.PDW==0);
        minP = 1e-300;
        LL = LL + sum(log(max(PrightHigh(dur_rightHigh),minP))) + sum(log(max(PrightLow(dur_rightLow),minP))) + ...
                  sum(log(max(PleftHigh(dur_leftHigh),minP)))   + sum(log(max(PleftLow(dur_leftLow),minP)));
        n = n + length(dur_rightHigh) + length(dur_rightLow) + length(dur_leftHigh) + length(dur_leftLow);
    end
    
    % copy trial params to data struct 'fit', then replace data with model vals
    fit = data;
    fit = rmfield(fit,'choice');
    try fit = rmfield(fit,'correct'); end %#ok<TRYNC>
    fit.pRight = expectedPright;
    if options.conftask==1 % SEP
        keyboard % not ready
    %     fit.conf = conf_model_trialwise;
    %     fit.PDW = nan(size(fit.conf));
    elseif options.conftask==2 % PDW
        fit.conf = expectedPhigh;
        fit.PDW = expectedPhigh; % legacy
        fit.pRightHigh = expectedPrightHigh;
        fit.pRightLow = expectedPrightLow;
    end
    if options.RTtask            
        fit.RT = expectedRT;
        fit.RThigh = expectedRThigh;
        fit.RTlow = expectedRThigh;
    end
    
    % also stored the 'parsed' values, for later plotting
    parsedFit = struct();
    parsedFit.pRight = pRight_model;
    if options.conftask==1 % SEP
        keyboard % not ready
%         parsedFit.confMean = meanConf_model;
    elseif options.conftask==2 % PDW
        parsedFit.pHigh = pHigh_model;        
        parsedFit.pHighCorr = pHighCorr_model;
        parsedFit.pHighErr = pHighErr_model;
        parsedFit.pRightHigh = pRightHigh_model;
        parsedFit.pRightLow = pRightLow_model;
    end
    if options.RTtask
        parsedFit.RTmean = meanRT_model;
        parsedFit.RTmeanHigh = meanRThigh_model;
        parsedFit.RTmeanLow = meanRTlow_model; 
    end

catch Error
    fprintf('\n\nERROR!!!\n'); fprintf('run %d\n', call_num);
    fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\talpha= %g\n', k, B, theta, alpha);
    fprintf('s:%d\n', s);
    for i = 1 : length(Error.stack)
        fprintf('entry %d of error stack\n', i);
        fprintf('\tfile: %s\n', Error.stack(i).file);
        fprintf('\tfunction: %s\n', Error.stack(i).name);
        fprintf('\tline: %d\n', Error.stack(i).line);
    end
    rethrow(Error);
end

err = -LL;
%     fit = data; % fit isn't needed here, but keep for consistency w 2D

%print progress report
if options.feedback
    fprintf('\n\n\n****************************************\n');
    fprintf('run %d\n', call_num);
    fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\talpha= %g\n', k, B, theta, alpha);
    fprintf('err: %f\n', err);
    fprintf('number of processed trials: %d\n', n);
end

if options.feedback==2 && strcmp(options.fitMethod,'fms')
    figure(options.fh);
    plot(call_num, err, '.','MarkerSize',14);
    drawnow;
end
call_num = call_num + 1;
toc

    
end



% retrieve the full parameter set given the adjustable and fixed parameters 
function param2 = getParam(param1, guess, fixed)
  param2(fixed==0) = param1(fixed==0);       %get adjustable parameters from param1
  param2(fixed==1) = guess(fixed==1);   %get fixed parameters from guess
end

