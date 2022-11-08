function [weights V lambda cp_model] = Weights(correlations,cp,use)

% function [weights V lambda cp_model] = Weights(correlations,cp,use)
%
% Weight reconstruction from input correlations and CPs for a large
% neuronal population (limit of infinitely many neurons)
%
% INPUTS
% correlations: correlation matrix
% cp: observed choice probabilities
% use: set of indices of eigenfunctions that are to be used in the
%      reconstruction (in descending order, i.e. use=[1 2] implies that 
%      only the eigenfunctions with the largest and second-largest
%      eigenvalue are to be used
%
% OUTPUTS
% weights: implied weights for a large neuronal population
% V: set of eigenfunctions used for the weight reconstruction
% lambda : eigenvalues of the eigenfunctions used for weight reconstruction
% cp_model: choice probabilities implied by the reconstructed weights

[V D]=eig(correlations);

[lambda idx]=sort(diag(D),'descend'); % sort by size of eigenvalue
V=V(:,idx); % apply the same sorting as for lambda's
nu=1./lambda(use).*(V(:,use)'*(cp-0.5));
weights=V(:,use)*nu;

if nargout>3
  % assuming constant variances across population for model
  cp_model=CP(correlations,weights,'infinite',ones(length(cp),1)); 
end