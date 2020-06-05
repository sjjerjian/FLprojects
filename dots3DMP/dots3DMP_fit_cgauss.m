%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

muPMF = nan(length(mods),length(cohs),length(deltas));
muPMFse = nan(length(mods),length(cohs),length(deltas));
sigmaPMF = nan(length(mods),length(cohs),length(deltas));
sigmaPMFse = nan(length(mods),length(cohs),length(deltas));

amplConf = nan(length(mods),length(cohs),length(deltas));
amplConfse = nan(length(mods),length(cohs),length(deltas));
muConf = nan(length(mods),length(cohs),length(deltas));
muConfse = nan(length(mods),length(cohs),length(deltas));
sigmaConf = nan(length(mods),length(cohs),length(deltas));
sigmaConfse = nan(length(mods),length(cohs),length(deltas));
baselineConf = nan(length(mods),length(cohs),length(deltas));
baselineConfse = nan(length(mods),length(cohs),length(deltas));

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
        beta = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
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
        beta = fminsearch(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,data.conf(I),data.heading(I)), [0.7 0 4 0.1]);
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
end