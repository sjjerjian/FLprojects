function y = logit(p,base)
% logit(p) returns log(p/(1-p))
% p can be a vector or a matrix
% useful for logistic regression.
% recall that p = exp(logit(p)) / (1+ exp(logit(p)))
% logit(p,b) uses b instead of e for the log base.
% So, a deciban is 10*logit(p,10)
%
% see logitinv for inverse. It also takes an optional second (base)
% argument.

% 1/29/98 mns wrote it
% 10/1/08 added 2nd arg


y = real(log(p ./ (1-p)));
if nargin>1
    y = y/log(base);
end

% I'm not sure the 'real' part was necessary, but there were situations
% encounterd in which matlab returned complex vals (the imag part was 0).  I
% have yet to understand why.
