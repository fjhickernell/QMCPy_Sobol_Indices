d = 2;
time_vec = (1:d)/d;
Sigma=bsxfun(@min, time_vec', time_vec);

% % Time differencing
fTD = @(x) max(sqrt(1/d)*cumsum(gail.stdnorminv(x'), 2))';
% Cholesky
A = chol(Sigma,'lower'); % A*A' = Sigma
fChol = @(x) max((A*(gail.stdnorminv(x)')))';
%PCA
[Eigenvectors,Eigenvalues]=eig(Sigma,'vector');
[~, order] = sort(Eigenvalues, 'descend');
A = Eigenvectors(:,order)*diag(Eigenvalues(order).^(1/2));% A*A' = Sigma
fPCA = @(x) max((A*(gail.stdnorminv(x)')))';

% Sobol indices estimation
hyperbox = [zeros(1,d) ; ones(1,d)];
abstol = 1e-3;
reltol = 0;
mmin = 9;
mmax = 22;
[q,app_int,out_param] = cubSobol_SI_all_g(fPCA,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax);
round(q*100)