function cp = CP(covariances,weights,FINITE,variances)

% function cp = CP(covariances,weights,FINITE,variances)
% computes choice probabilities cp implied by covariances and weights
%
% if FINITE='finite' (default) a finite population is assumed
%
% if FINITE='infinite' an infinitely large population is assumed where
% covariances and weights contain the discretization of the covariance
% function and the weights function
% 
% currently works only for 2D C and 1D weights but can be generalized to
% incorporate time or, in the case of infinitely large populations,
% multiple dependencies

if nargin<3, FINITE='finite'; end

weights=weights(:); % make sure it's a vector
switch FINITE
  case 'finite', variances=diag(covariances);
  case 'infinite', variances=variances(:);
  otherwise
    error(FINITE);
end

xsi=(covariances*weights)./sqrt(variances*(weights'*covariances*weights));
cp=0.5+2/pi*sign(xsi).*atan(1./sqrt(2./xsi.^2-1)); % equation (1) in paper
