function [error] = logitFitting(independent,dependent, betas)


% calculate probability using a model (Logit for now)
p = 1./(1+exp(-(betas' * independent'))); %should give you many values for P, a P for each trial
%p(p == 0) = eps;
p = p - (p > (1.0-eps))*eps + (p < eps)*eps; %analyze what this does, this solves any thing that goes to NaN or Inf
%Then calculate the loglikehood and min this
%keyboard
logll = sum(log(p).*(dependent' == 1))  + sum(log(1-p).*(dependent' == 0));

error = -logll;
end

