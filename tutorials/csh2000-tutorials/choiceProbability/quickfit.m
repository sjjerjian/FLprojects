function [alpha, beta, llik,abse,theta0] = quickfit(data)
% QUICKFIT fits a weibull function to data using maximum likelihood
%	maximization under binomial assumptions.  It uses FITFUNW for
%	error calculation.  Data must be in 3 columns: x, %-correct, Nobs
%	usage: [alpha, beta,llik,abse] = quickfit(data)
%		alpha and beta are the threshold and slope parameters. llik is the
%		log likelihood of obtaining data given the fit. abse is a
%		2-vector with standard errors for alpha and beta

% update 4/1/95 to use optimization toolbox, contrained fit -- nope
%  old code is quickfit_no_opt.m
% 12/11/96  jdr & mns fixed guessing bug (we hope) 
% 10/16/02 mns updated initial guess using a crude logistic regression.
% 3/21/03 mns fixed bug in logistic guess and added SE 
% 6/24/09 mns bypassed fminunc (needs real update)

% generate threshold guess using logisitic (update 3/21/03)
q = ones(2,1);
pobs = data(:,2);
% Get rid of 1's and 0's
pobs(pobs > 1-eps) = nan;
pmax = nanmax(pobs) + .9* (1 - nanmax(pobs));
pobs(isnan(pobs)) = pmax;
pobs(pobs < eps) = nan;
pmin = .1 * nanmin(pobs);
pobs(isnan(pobs)) = pmin;
% convert to logit
logLR = logit(pobs);
% do least squares regression
blogistic = regress(logLR,[ones(size(data,1),1) data(:,1)]);
% What p correspondes to x=alpha?
pAtAlpha = 1 - 0.5*exp(-1);
% take the inverse to get alpha
q(1,1) = (logit(pAtAlpha)-blogistic(1)) / blogistic(2);

% Use nelder-mead 1st
[theta0,fval0] = fminsearch(@(q) fitfunw(q,data),q);

if nargout > 3
    % refine with fminunc to get hessian
    [theta,fval,xflg,outp,grd,hess] = fminunc(@(q) fitfunw(q,data),theta0);

    alpha = theta(1,1);
    beta = theta(2,1);
    llik = fval;
    abse = sqrt(diag(inv(hess)));
else
    alpha = theta0(1);
    beta = theta0(2);
    llik = fval0;
end


