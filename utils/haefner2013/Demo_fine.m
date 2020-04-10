clear all;
figure;
% n_neurons  =2000; % number of neurons in population
n_neurons  =46; % number of neurons in population
obs_reconst=1:46; % observed neurons
% obs_reconst=5:20:2000; % observed neurons
disp(['using ' num2str(numel(obs_reconst)) ' neurons for reconstruction']);
n_reconst  =100;  % number of discretization points for reconstruction
use_reconst=1:5;  % use eigenvectors with 5 largest eigenvalues for reconstruction

% Define limited range correlations between 0.0 and 0.3:
% c=linspace(0.3,0.0,n_neurons/2); c=[c fliplr(c)];
% Cor=toeplitz(c); % correlation function

Cor = make_cormat;

% Cor = Cor-0.125;
% figure; contourf(Cor); colorbar;


% Assume heterogenous population with mean response variance of 10 and
% minimum variance of 1
v=1+poissrnd(9,[n_neurons 1]);
Cov=sqrt(diag(v))*Cor*sqrt(diag(v)); % Covariances from correlations

subplot(2,2,1);  hold on; title('weights');
% define uniform weights in each pool
w=[ones(1,n_neurons/2) -ones(1,n_neurons/2)]'; 
plot(w/mean(w(1:n_neurons/2)),'b-'); % plot (normalized) weights

subplot(2,2,2);  hold on; title('choice probability');
cp=CP(Cov,w,'finite'); % computation of actual choice probabilities
plot(cp,'b-');

% use only a subset of all neurons for reconstruction
cp=cp(obs_reconst);
Cov=Cov(obs_reconst,obs_reconst);
% adding noise to choice probabilities simulating 200 trials
n_trials=200;
cp=binornd(n_trials,cp)/n_trials;
subplot(2,2,2); hold on;
scatter(obs_reconst,cp,5,'r');
legend('ground truth','observation');
% reconstruction of weights based on eigenvectors with 5 largest eigenvalues
Cor=corrcov(Cov); % compute correlations from covariances for reconstruction
[weights,V,lambda,cp_model]=Weights(Cor,cp,use_reconst); 

subplot(2,2,1);
plot(obs_reconst,weights/mean(weights(1:length(obs_reconst)/2)),'r-');
legend('ground truth','reconstruction');
subplot(2,2,3); hold on; title('eigenvalues');
plot(lambda,'b.');
plot(lambda(use_reconst),'r.');
legend('all EVs','used for reconstruction');
subplot(2,2,4); hold on; title('eigenvectors used for reconstruction');
plot(V(:,use_reconst));
for i=1:numel(use_reconst)
  leg{i}=['EV=' num2str(lambda(use_reconst(i)),2)];
end
legend(leg);
