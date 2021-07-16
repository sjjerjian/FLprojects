function [betas, error] = binarylogll_logit(independent, dependent)

betas0 = rand(size(independent, 2) + 1, 1);
independent = [ones(size(independent,1),1) independent]; %Add the bias term

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun' , 1e-7, 'TolX', 1e-7);
betas = fminsearch(@(betas0) logitFitting(independent, dependent, betas0), betas0, options);

error = logitFitting(independent, dependent, betas);
end

