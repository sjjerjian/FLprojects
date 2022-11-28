function gfit = dots3DMP_fit_cgauss_NN(data,mods,cohs,deltas)

% define anonymous functions for fitting:

cgauss = @(b,hdg) 1/2 * ( 1 + erf( (hdg-b(1))./(b(2)*sqrt(2)) ) );
    % for probabilities, error is negative log likelihood of observing the data, which is
    % [ log(Pright(hdg)) + log(1-(~Pright(hdg))) ]
cgauss_err = @(param,choice,hdg) -(sum(log(cgauss(param,hdg(choice))))+sum(log(1-cgauss(param,hdg(~choice))))); 

unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)

%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

% initialize vars for storing param fits
n = nan(length(mods),length(cohs),length(deltas));

muPMF = n; muPMFse = n;
sigmaPMF = n; sigmaPMFse = n;

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
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I)==1,data.heading(I)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==1,data.heading(I)), [0 3]);
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
    
end


%% now separate by delta

for c = 1:length(cohs)
    % choice
    for d = 1:length(deltas)     % m c d h
        I = data.modality==3 & data.coherence==cohs(c) & data.delta==deltas(d);
        beta = fminsearch(@(x) cgauss_err(x,data.choice(I)==1,data.heading(I)), [0 3]);
        [betaUnc,~,~,~,~,hessian] = fminunc(@(x) cgauss_err(x,data.choice(I)==1,data.heading(I)), [0 3]);
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
    
end

gfit = struct();
gfit.muPMF = muPMF;
gfit.muPMFse = muPMFse;
gfit.sigmaPMF = sigmaPMF;
gfit.sigmaPMFse = sigmaPMFse;
gfit.func.cgauss = cgauss;
gfit.func.cgauss_err = cgauss_err;
