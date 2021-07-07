function [tuning_curves, posterior, pop_lh, xq, simChoice] = dots3DMP_likelihood_decoding(data,meanFRs,hdgs,mods,cohs,deltas,numtrs,step)

% numtrs - number of simulated trials
% step - tuning curve stepsize

lh_func = @(tc,r_obs) (exp(-tc) .* bsxfun(@power,tc,r_obs)) ./ gamma(r_obs); % vectorize it

numunits = length(unique(data.unitnum));
xq = min(hdgs):step:max(hdgs);

% pre-allocate
tuning_curves = nan(length(mods),length(cohs),length(deltas),length(xq),numunits);
r_obs         = nan(length(mods),length(cohs),length(deltas),length(hdgs),numunits,numtrs);
lh            = nan(length(mods),length(cohs),length(deltas),length(hdgs),length(xq),numunits,numtrs);

% surely there's a way to vectorize this...anyway, loops for now
for m=1:length(mods)
for c=1:length(cohs)
for d=1:length(deltas)
for u=1:numunits
    x = squeeze(meanFRs(m,c,d,:,u));
        
    if any(isnan(x)),continue,end % skip

    % linear interpolation
    tuning_curves(m,c,d,:,u) = interp1(hdgs,x,xq);
    
    % linear regression
%     p = polyfit(hdgs,x,1);
%     y = polyval(p,xq);
%     tuning_curves(m,c,d,:,u) = y;


    % random draws of single-trial firing for each condition
    % ... and compute likelihood 
    
    % tuning curve for given stimulus
    tc = squeeze(tuning_curves(m,c,deltas==0,:,u)); % use delta==0
    
    for h=1:length(hdgs)
        
        % poisson random draw from mean interpolated response
%         [~,idx] = min(abs(xq-hdgs(h)));
%         r = poissrnd(tuning_curves(m,c,d,idx,u),numtrs,1);  

        % can't we just use the meanFRs? 
        r = poissrnd(meanFRs(m,c,d,h,u),numtrs,1);
        r_obs(m,c,d,h,u,:) = r;
        
        lh(m,c,d,h,:,u,:) = lh_func(tc,r');
        
    end
end
end
end
end

lh(lh==0) = 1e-3;               % kluge to avoid prod(0...)
pop_lh = squeeze(prod(lh,6));   % population likelihood, product across units

% normalize to posterior, assuming flat prior
posterior = pop_lh ./ sum(pop_lh,5);

%%
[~,zerohdg] = min(abs(xq));
negSum = squeeze(nansum(posterior(:,:,:,:,1:zerohdg-1,:,:),5));
posSum = squeeze(nansum(posterior(:,:,:,:,zerohdg+1:end,:,:),5));

simChoice = posSum > negSum;
