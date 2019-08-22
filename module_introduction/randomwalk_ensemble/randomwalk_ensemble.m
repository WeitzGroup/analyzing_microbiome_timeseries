% Generate an ensemble of random walks and compute correlations between
% pairs. Results from this data are shown in figure_autoregression, 3rd
% column.

clear;

fprintf('generating ensemble of random walks...\n');

% setup
N = 20000; % number of timeseries (# of correlations = N/2)
NT = 100; % number of timepoints in each timeseries

% generate random walks
epsilon = 2*rand(NT,N)-1;
X = cumsum(epsilon);

fprintf('computing residuals...\n');

% compute residual timeseries
X1 = X(1:end-1,:);
X2 = X(2:end,:);
phi1_X = sum((X1-mean(X1)).*(X2-mean(X2)))./sum(X1.^2-mean(X1).^2);
Xresidual = X2 - X1.*phi1_X;

fprintf('computing correlations...\n');

% no cross-correlations (only use a time-series once)
rho_X = zeros(N/2,1);
pval_X = zeros(N/2,1);
rho_Xresidual = zeros(N/2,1);
pval_Xresidual = zeros(N/2,1);
for i = 1:N/2
    ID1 = i*2-1;
    ID2 = i*2;
	[rho_X(i),pval_X(i)] = corrcoef(X(:,ID1),X(:,ID2)); 
    [rho_Xresidual(i),pval_Xresidual(i)] = corr(Xresidual(:,ID1),Xresidual(:,ID2));
end

fprintf('saving results...\n');

rho = rho_X;
save('randomwalk_ensemble_originals','rho')

rho = rho_Xresidual;
save('randomwalk_ensemble_residuals','rho');

fprintf('randomwalk_ensemble complete\n');

clear;