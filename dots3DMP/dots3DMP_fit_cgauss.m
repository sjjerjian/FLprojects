% define anonymous functions for fitting:

cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
    % for probabilities, error is negative log likelihood of observing the data, which is
    % [ log(Pright(hdg)) + log(1-(~Pright(hdg))) ]
cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice))))); 

% fit a (flipped) gaussian for conf and RT for now, until DDM fits are good
% b(2) and b(3) are mu and sigma
% b(1) is 'amplitude' above baseline, b(4) is baseline level

% 'baseline' for flippedGauss is the highBet side, because of flip
if conftask==1 % continuous, sacc endpoint    
    % for continuous values, error is sum squared error
    flippedGauss = @(b,hdg) 1 - ( b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4));
    flippedGauss_err = @(param,SEP,hdg) nansum((flippedGauss(param,hdg)-SEP).^2);
elseif conftask==2 % PDW, probabilities
    flippedGauss = @(b,hdg) 1 - ( min(max(b(1),0),1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + min(max(b(4),0),1));

    % negative log likelihood of observing PDW data
    % log prob of observing high bet on all trials where subj bet
    % high, + log prob of low bet on trials where subj bet low
    flippedGauss_err = @(param,pdw,hdg) -( sum(log(flippedGauss(param,hdg(pdw)))) + sum(log(1-flippedGauss(param,hdg(~pdw)))) );
end

% RT
gauss = @(b,hdg) b(1) .* exp(-(hdg-b(2)).^2 ./ (2*b(3).^2)) + b(4);
gauss_err = @(param,SEP,hdg) sum((gauss(param,hdg)-SEP).^2);

unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)

% parameter initial guesses
guess_fgauss = [0.2 0 6 0];
guess_gauss = [1 0 6 1];


%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

% initialize vars for storing param fits
% deal func looks nicer, but is slow for some reason...
n = nan(length(mods),length(cohs),length(deltas));

muPMF = n; muPMFse = n;
sigmaPMF = n; sigmaPMFse = n;

amplConf = n; amplConfse = n;
muConf = n; muConfse = n;
sigmaConf = n; sigmaConfse = n;
baselineConf = n; baselineConfse = n;

amplRT = n; amplRTse = n;
muRT = n; muRTse = n;
sigmaRT = n; sigmaRTse = n;
baselineRT = n; baselineRTse = n;

for c = 1:length(cohs)
    % choice
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
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(m,c,D) = SE(1);
        sigmaPMFse(m,c,D) = SE(2);
        if unc
            muPMF(m,c,D) = betaUnc(1);
            sigmaPMF(m,c,D) = betaUnc(2);
        else
            muPMF(m,c,D) = beta(1);
            sigmaPMF(m,c,D) = beta(2);
        end
    end

    % conf
    for m = 1:length(mods)        
        if m==1
            I = data.modality==mods(m);
        else
            if D==length(deltas)+1
                I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
            else
                I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
            end
        end
        
        if conftask==1 % sacc endpoint
            
            beta = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss);
        elseif conftask==2 % PDW
            beta = fminsearch(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss);
        end
            
        SE = sqrt(diag(inv(hessian)));
        amplConfse(m,c,D) = SE(1);
        muConfse(m,c,D) = SE(2);
        sigmaConfse(m,c,D) = SE(3);
        baselineConfse(m,c,D) = SE(4);
        if unc
            amplConf(m,c,D) = betaUnc(1);
            muConf(m,c,D) = betaUnc(2);
            sigmaConf(m,c,D) = betaUnc(3);
            baselineConf(m,c,D) = betaUnc(4);
        else
            amplConf(m,c,D) = beta(1);
            muConf(m,c,D) = beta(2);
            sigmaConf(m,c,D) = beta(3);
            baselineConf(m,c,D) = beta(4);
        end
    end
    
    % RT
    if ~isnan(RTmean(1,1,2))
        for m = 1:length(mods)        
            if m==1
                I = data.modality==mods(m);
            else
                if D==length(deltas)+1
                    I = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
                else
                    I = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(D);
                end
            end
            beta = fminsearch(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(m,c,D) = SE(1);
            muRTse(m,c,D) = SE(2);
            sigmaRTse(m,c,D) = SE(3);
            baselineRTse(m,c,D) = SE(4);
            if unc
                amplRT(m,c,D) = betaUnc(1);
                muRT(m,c,D) = betaUnc(2);
                sigmaRT(m,c,D) = betaUnc(3);
                baselineRT(m,c,D) = betaUnc(4);
            else
                amplRT(m,c,D) = beta(1);
                muRT(m,c,D) = beta(2);
                sigmaRT(m,c,D) = beta(3);
                baselineRT(m,c,D) = beta(4);
            end
        end
    end
    
end


%% now separate by delta

for c = 1:length(cohs)
    % choice
    for d = 1:length(deltas)     % m c d h
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==2,data.heading(I)), [0 3]);
        SE = sqrt(diag(inv(hessian)));
        muPMFse(3,c,d) = SE(1);
        sigmaPMFse(3,c,d) = SE(2);        
        if unc
            muPMF(3,c,d) = betaUnc(1);
            sigmaPMF(3,c,d) = betaUnc(2);
        else
            muPMF(3,c,d) = beta(1);
            sigmaPMF(3,c,d) = beta(2);
        end
    end

    % conf
    for d = 1:length(deltas)
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        
        if conftask==1 % sacc endpoint
            beta = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), guess_fgauss);
        elseif conftask==2 % PDW
            beta = fminsearch(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.PDW(I)==1,data.heading(I)), guess_fgauss);
        end
            
        SE = sqrt(diag(inv(hessian)));
        amplConfse(3,c,d) = SE(1);
        muConfse(3,c,d) = SE(2);
        sigmaConfse(3,c,d) = SE(3);
        baselineConfse(3,c,d) = SE(4);
       if unc
            amplConf(3,c,d) = betaUnc(1);
            muConf(3,c,d) = betaUnc(2);
            sigmaConf(3,c,d) = betaUnc(3);
            baselineConf(3,c,d) = betaUnc(4);
        else
            amplConf(3,c,d) = beta(1);
            muConf(3,c,d) = beta(2);
            sigmaConf(3,c,d) = beta(3);
            baselineConf(3,c,d) = beta(4);
        end
    end
    
    % RT
    if ~isnan(RTmean(1,1,2))
        for d = 1:length(deltas)
            I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
            beta = fminsearch(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss);
            [betaUnc,~,~,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I),data.heading(I)), guess_gauss);
            SE = sqrt(diag(inv(hessian)));
            amplRTse(3,c,d) = SE(1);
            muRTse(3,c,d) = SE(2);
            sigmaRTse(3,c,d) = SE(3);
            baselineRTse(3,c,d) = SE(4);
           if unc
                amplRT(3,c,d) = betaUnc(1);
                muRT(3,c,d) = betaUnc(2);
                sigmaRT(3,c,d) = betaUnc(3);
                baselineRT(3,c,d) = betaUnc(4);
            else
                amplRT(3,c,d) = beta(1);
                muRT(3,c,d) = beta(2);
                sigmaRT(3,c,d) = beta(3);
                baselineRT(3,c,d) = beta(4);
            end
        end
    end
    
end