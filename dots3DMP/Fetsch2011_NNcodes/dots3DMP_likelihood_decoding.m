function [tuning_curves, posterior, pop_lh, hdg_intp, simChoice] = dots3DMP_likelihood_decoding(data,meanFRs,hdgs,mods,cohs,deltas,numtrs,step)

% numtrs - number of simulated trials
% step - tuning curve stepsize

lh_func = @(g,tc,r_obs) (exp(-g*tc) .* bsxfun(@power,g*tc,r_obs)) ./ factorial(r_obs); % vectorized

numunits = length(unique(data.unitnum));
hdg_intp = min(hdgs):step:max(hdgs);

% pre-allocate
tuning_curves = nan(length(mods),length(cohs),length(deltas),length(hdg_intp),numunits);
r_obs         = nan(length(mods),length(cohs),length(deltas),length(hdgs),numunits,numtrs);
lh            = nan(length(mods),length(cohs),length(deltas),length(hdgs),length(hdg_intp),numunits,numtrs);

g = 1; % Poisson gain

% surely there's a way to vectorize this...anyway, loops for now
for m=1:length(mods)
for c=1:length(cohs)
for d=1:length(deltas)
for u=1:numunits
    x = squeeze(meanFRs(m,c,d,:,u));
        
    if any(isnan(x)),continue,end % skip

    % linear interpolation
    tuning_curves(m,c,d,:,u) = interp1(hdgs,x,hdg_intp,'linear');
    
    % linear regression
%     p = polyfit(hdgs,x,1);
%     y = polyval(p,xq);
%     tuning_curves(m,c,d,:,u) = y;

    % random draws of single-trial firing for each condition
    % ... and compute likelihood 
    
    % tuning curve for given stimulus
    tc = squeeze(tuning_curves(m,c,deltas==0,:,u)); % use delta==0!!
    
    for h=1:length(hdgs)
        
        J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d) & data.unitnum==u;
        
        r_trials = round(data.spikeRates(J));
        if isempty(r_trials), continue, end
        shuf_ind = randi(length(r_trials),1,numtrs);

        r = r_trials(shuf_ind);
        r_obs(m,c,d,h,u,:) = r;
        
        if any(sum(bsxfun(@power,g*tc,r')) == Inf) || any(factorial(r) == Inf)
            %                 disp(['Inf problem! Cell ID #' num2str(n)]);
            tc = tc .* 0.8;  % This kluges one cell whose high FR makes the equation exceed
            r = round(r.*0.8); % the largest representable positive floating point number
        end
        if any(sum(bsxfun(@power,g*tc,r')) == Inf) || any(factorial(r) == Inf)
            %             disp(['Inf problem after scale-down! Cell ID #' num2str(n)]);
            tc = tc .* 0.5;  % One(?) cell needs even more of a scale-down
            r = round(r.*0.5);
        end
        if any(sum(bsxfun(@power,g*tc,r')) == Inf) || any(factorial(r) == Inf)
            disp(['Inf problem after double scale-down! Cell ID #' num2str(n)]);
            keyboard
        end
        if sum(tc) == 0
            disp(['Zero problem: Cell ID #' num2str(n)]);
            tc = ones(size(tc)); % cell does not contribute, but does not zero out the likelihood either
        end
        
        lh(m,c,d,h,:,u,:) = lh_func(g,tc,r');
        
    end
end
end
end
end

lh(lh==0) = .01;               % kluge to avoid prod(0...)
pop_lh = squeeze(prod(lh,6));   % population likelihood, product across units

% normalize to get posterior, assuming flat prior
posterior = pop_lh ./ sum(pop_lh,5);

%%

% [~,zerohdg] = min(abs(hdg_intp));
zerohdg = find(hdg_intp==0);
negSum = squeeze(nansum(posterior(:,:,:,:,1:zerohdg-1,:,:),5));
posSum = squeeze(nansum(posterior(:,:,:,:,zerohdg+1:end,:,:),5));
simChoice = posSum > negSum * 2 - 1;

% [~,maxpos] = max(posterior,[],5);
% MAPhdg     = squeeze(hdg_intp(maxpos));
% simChoice  = sign(MAPhdg);

simChoice(simChoice==0) = sign(randn);

