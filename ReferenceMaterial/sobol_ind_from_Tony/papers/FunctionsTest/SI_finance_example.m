clearvars

abstol = 1e-3;
reltol = 0; % Pure absolute tolerance
mmax = 22; % I adjust that not to run out of memory. It can go up to 54. Type help cubSobol_SI_g for more information.


% %% European Call with PCA
% S0 = 100; K = 100; sgma = 0.05; r = sgma^2/2; T = 1; price = normcdf((log(S0/K)+(r+sgma^2/2)*T)/(sgma*sqrt(T)))*S0 - normcdf((log(S0/K)+(r-sgma^2/2)*T)/(sgma*sqrt(T)))*K*exp(-r*T);
% 
% d = 8;
% time_disc = linspace(T/d,T,d);
% Sigma=bsxfun(@min,time_disc',time_disc);  % Sigma(i,j) = min(time_disc(i),time_disc(j))
% [Evec,Eval]=eig(Sigma,'vector'); % Jordan decomposition of Sigma
% [~,order] = sort(Eval,'descend');
% A = Evec(:,order)*diag(Eval(order).^(1/2)); % A*A' = Sigma. PCA substitution
% 
% 
% f = @(x) max(0,(S0*exp((r-sgma^2/2)*bsxfun(@plus,time_disc,zeros(size(x,1),1))+sgma*norminv(x)*A'))*[zeros(d-1,1);1]-K);
% hyperbox = [zeros(1,d) ; ones(1,d)];
% 
% j = 1;
% [SI, approx_price, out_param] = cubSobol_SI_g(f,hyperbox,j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
% SI
% approx_price

%% Geometric Asian Option with PCA
S0 = 100; K = 100; sgma = 0.05; r = 0.02; T = 2; price = 000;

d = 6;
time_disc = linspace(T/d,T,d);
Sigma=bsxfun(@min,time_disc',time_disc);  % Sigma(i,j) = min(time_disc(i),time_disc(j))
[Evec,Eval]=eig(Sigma,'vector'); % Jordan decomposition of Sigma
[~,order] = sort(Eval,'descend');
A = Evec(:,order)*diag(Eval(order).^(1/2)); % A*A' = Sigma. PCA substitution

f = @(x) max(0,S0*prod(exp((r-sgma^2/2)*bsxfun(@plus,time_disc,zeros(size(x,1),1))+sgma*norminv(x)*A'),2).^(1/d)-K);
hyperbox = [zeros(1,d) ; ones(1,d)];

tic
for j=1:d
    [si, approx_price, out_param] = cubSobol_SI_g(f,hyperbox,j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    SI(1,j) = si;
    [si, approx_price, out_param] = cubSobol_SI_g(f,hyperbox,-j,'abstol',abstol,'reltol',reltol,'mmax',mmax);
    SI(2,j) = si;
end
toc
approx_price