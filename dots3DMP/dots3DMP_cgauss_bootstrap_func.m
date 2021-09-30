function [gfitBoot,wvesBoot] = ...
    dots3DMP_cgauss_bootstrap_func(data,gfit,mods,cohs,deltas,nboots)
% runs exactly the same gaussian fit as the regular cgauss, and computes
% weights like the weight

% if gfit passed in, use those funcs, otherwise redefine them!
if ~isempty(gfit)
    keyboard
    cgauss = gfit.choice.func;
    cgauss_err = gfit.choice.err;
    guess_cgauss = gfit.choice.guess; 
    
    flippedGauss = gfit.conf.func;
    flippedGauss_err = gfit.conf.err;
    guess_fgauss = gfit.conf.guess; 
    
    gauss = gfit.RT.func;
    cgauss_err = gfit.RT.err;
    guess_gauss = gfit.RT.guess; 
else
    cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
    cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice)))));
    
    % CONFIDENCE - 'flipped' Gaussian
    if conftask==1 % continuous, sacc endpoint
        flippedGauss = @(b,hdg) 1 - ( b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));
        flippedGauss_err = @(param,SEP,hdg) nansum((flippedGauss(param,hdg)-SEP).^2);
    elseif conftask==2 % PDW, probabilities
        flippedGauss = @(b,hdg) 1 - ( b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));
        flippedGauss_err = @(param,pdw,hdg) -( sum(log(flippedGauss(param,hdg(pdw)))) + sum(log(1-flippedGauss(param,hdg(~pdw)))) );
    end
    
    % RT - Gaussian, error is sum squared because RT is cont variable
    gauss = @(b,hdg) b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4);
    gauss_err = @(param,RT,hdg) sum((gauss(param,hdg)-RT).^2);
        
    % parameter initial guesses
    guess_cgauss = [0 3];
    guess_fgauss = [0.2 0 6 0.5];
    guess_gauss  = [0.5 0 6 1];
end

if nboots==0, return, end

gfitBoot.choice.mu      = cell(nboots,1);
gfitBoot.choice.sigma   = cell(nboots,1);
gfitBoot.conf.ampl      = cell(nboots,1);
gfitBoot.conf.mu        = cell(nboots,1);
gfitBoot.conf.sigma     = cell(nboots,1);

unc = 0;
options = optimset('Display','off');

Nmcd = nan(length(mods),length(cohs),length(deltas));

if conftask==1
    conf = data.conf;
elseif conftask==2
    conf = data.PDW;
end

for n = 1:nboots
    
    if mod(n,10)==0
        fprintf('\nprogress: %d out of %d\n', n, nboots);
    end
    
    gfitBoot.choice.mu{n}    = Nmcd;
    gfitBoot.choice.sigma{n} = Nmcd;
    gfitBoot.conf.ampl{n}    = Nmcd;
    gfitBoot.conf.mu{n}      = Nmcd;
    gfitBoot.conf.sigma{n}   = Nmcd;
    
    % first the single-cues, and comb with either:
    % all trials irrespective of delta
    D = length(deltas)+1; % (the extra column we made for pooling across deltas)
    % OR select just delta=0:
    % D = find(deltas==0);
    
    for c = 1:length(cohs)
        for m = 1:length(mods)     % m c d h
            if m==1
                I = data.modality==mods(m);
            else
                if D==length(deltas)+1
                    I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
                else
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
                end
            end
            
            bootI = randsample(find(I),sum(I),true);
            
            if unc
                [beta,~,~,~,~,~] = fminunc(@(x) cgauss_err(x,data.choice(bootI)==2,data.heading(bootI)), [0 3], options);
                if conftask>0
                    [betaConf,~,~,~,~,~] = fminunc(@(x) flippedGauss_err(x,conf(bootI),data.heading(bootI)), [0.7 0 4 0.1], options);
                end
                
            else
                beta = fminsearch(@(x) cgauss_err(x,data.choice(bootI)==2,data.heading(bootI)), [0 3], options);
                if conftask>0
                    betaConf = fminsearch(@(x) flippedGauss_err(x,conf(bootI),data.heading(bootI)), [0.7 0 4 0.1], options);
                end
            end
            gfitBoot.choice.mu{n}(m,c,D)     = beta(1);
            gfitBoot.choice.sigma{n}(m,c,D)  = beta(2);
            gfitBoot.conf.ampl{n}(m,c,D)  = betaConf(1);
            gfitBoot.conf.mu{n}(m,c,D)    = betaConf(2);
            gfitBoot.conf.sigma{n}(m,c,D) = betaConf(3);
        end
    end
    
    % then the comb only, separated by delta, to get biases
    for c = 1:length(cohs)
        for d = 1:length(deltas)     % m c d h
            I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
            
            bootI = randsample(find(I),sum(I),true);
            
            if unc
                [beta,~,~,~,~,~] = fminunc(@(x) cgauss_err(x,data.choice(bootI)==2,data.heading(bootI)), [0 3], options);
                if conftask>0
                    [betaConf,~,~,~,~,~] = fminunc(@(x) flippedGauss_err(x,conf(bootI),data.heading(bootI)), [0.7 0 4 0.1], options);
                end
            else
                beta = fminsearch(@(x) cgauss_err(x,data.choice(bootI)==2,data.heading(bootI)), [0 3], options);
                if conftask>0
                    betaConf = fminsearch(@(x) flippedGauss_err(x,conf(bootI),data.heading(bootI)), [0.7 0 4 0.1], options);
            
                end
            end
            
            gfitBoot.choice.mu{n}(3,c,d)     = beta(1);
            gfitBoot.choice.sigma{n}(3,c,d)  = beta(2);
            gfitBoot.conf.ampl{n}(3,c,d)  = betaConf(1);
            gfitBoot.conf.mu{n}(3,c,d)    = betaConf(2);
            gfitBoot.conf.sigma{n}(3,c,d) = betaConf(3);
        end
    end
    
    % lastly compute the weights
    for c = 1:length(cohs)      % m c d
        wvesBoot.choice.pred(n,c) = (1/gfitBoot.choice.sigma{n}(1,1,D)^2) / ((1/gfitBoot.choice.sigma{n}(1,1,D)^2) + (1/sigmaPMFboot{n}(2,c,D)^2));
        % m c d
        actual(1) = (gfitBoot.choice.mu{n}(3,c,1)-gfitBoot.choice.mu{n}(3,c,2)+(deltas(1)/2)) / deltas(1);
        actual(2) = (gfitBoot.choice.mu{n}(3,c,3)-gfitBoot.choice.mu{n}(3,c,2)+(deltas(3)/2)) / deltas(3);
        wvesBoot.choice.emp(n,c) = mean(actual);
        
        actual(1) = (gfitBoot.conf.mu{n}(3,c,1)-gfitBoot.conf.mu{n}(3,c,2)+(deltas(1)/2)) / deltas(1);
        actual(2) = (gfitBoot.conf.mu{n}(3,c,3)-gfitBoot.conf.mu{n}(3,c,2)+(deltas(3)/2)) / deltas(3);
        wvesBoot.choice.emp(n,c) = mean(actual);
    end
    
    
end