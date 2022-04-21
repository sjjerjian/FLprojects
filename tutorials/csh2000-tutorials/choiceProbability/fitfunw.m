function err = fitfunw(q,data)
%FITFUNW   Used by QUICKFIT.
%	FITFUNW(obj,q,data) returns the error between the data and the
%	values computed by the current function of weibul params.
%	FITFUNW assumes a function of the form
%
%	  y = 1 - .5 * exp( -(x/q(1))^q(2) )
%
%	where q(1) is alpha, and q(2) is beta.
%
%	The data is in columns such that 
%   data(:,1) is abscissa,
%	data(:,2) is observed percent correct (0..1), and
%	data(:,3) is number of observations.
%
%	The value of err is the -log likelihood of obtaining data
%	given the parameters q.
%

x = data(:,1);
y = data(:,2);
n = data(:,3);
z = 1 - .5 * exp( -(x/q(1)).^q(2) );
z = z - (z > (1.0-eps))*eps + (z < eps)*eps; 
llik = n .* y .* log(z) +  n .* (1-y) .* log(1-z);

% for moment, just return the square error
% err = norm(z-y);
err = -sum(llik);

