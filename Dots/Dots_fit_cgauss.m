% define anonymous functions for fitting:

cgauss = @(b,coh) 1/2 * ( 1 + erf( (coh-b(1))./(b(2)*sqrt(2)) ) );
    % for probabilities, error is negative log likelihood of observing the data, which is
    % [ log(Pright(coh)) + log(1-(~Pright(coh))) ]
cgauss_err = @(param,choice,coh) -(sum(log(cgauss(param,coh(choice))))+sum(log(1-cgauss(param,coh(~choice))))); 


% fit a (flipped) gaussian for PDW and RT for now, until DDM fits are good
flippedGauss = @(b,coh) 1 - ( min(max(b(1),0),1) .* exp(-(coh-b(2)).^2 ./ (2*b(3).^2)) + b(4));
flippedGauss_err = @(param,PDW,coh) -(sum(log(flippedGauss(param,coh(PDW))))+sum(log(1-flippedGauss(param,coh(~PDW))))); 

% % %quantify each Psure function (stimulated and unstimulated) with a bell-shaped function
% % psure_fcn_b3 = @(b,coh) min(max(b(1),0),1).*exp(-(coh+b(3)).^2 ./ (2*b(2).^2));
% % err_fcn_b3 = @(param,choice,coh) -(sum(log(psure_fcn_b3(param,coh(choice))))+sum(log(1-psure_fcn_b3(param,coh(~choice)))));


gauss = @(b,coh) b(1) .* exp(-(coh-b(2)).^2 ./ (2*b(3).^2)) + b(4);
    % for continuous values, error is sum squared error
gauss_err = @(param,RT,coh) sum((gauss(param,coh)-RT).^2);


unc = 0; % saves biases from fminunc instead of fminsearch (SEs always are fminunc, and plots are always fminsearch)

% parameter intial guesses
%              [A  mu sigma offset]
guess_fgauss = [0.2 0 0.1 0.1]; % PDW
guess_gauss =  [0.4 0 0.3 0.3]; % RT

%%

% conf
I = ~isnan(data.PDW);
beta = fminsearch(@(x) flippedGauss_err(x,logical(data.PDW(I)),data.coherence(I)), guess_fgauss);
[betaUnc,~,~,~,~,hessian] = fminunc(@(x) flippedGauss_err(x,logical(data.PDW(I)),data.coherence(I)), guess_fgauss);
SE = sqrt(diag(inv(hessian)));
amplConfse = SE(1);
muConfse = SE(2);
sigmaConfse = SE(3);
baselineConfse = SE(4);
amplConf = beta(1);
muConf = beta(2);
sigmaConf = beta(3);
baselineConf = beta(4);

% RT
I = ~isnan(data.RT);
beta = fminsearch(@(x) gauss_err(x,data.RT(I),data.coherence(I)), guess_gauss);
[betaUnc,~,~,~,~,hessian] = fminunc(@(x) gauss_err(x,data.RT(I),data.coherence(I)), guess_gauss);
SE = sqrt(diag(inv(hessian)));
amplRTse = SE(1);
muRTse = SE(2);
sigmaRTse = SE(3);
baselineRTse = SE(4);
amplRT = beta(1);
muRT = beta(2);
sigmaRT = beta(3);
baselineRT = beta(4);

