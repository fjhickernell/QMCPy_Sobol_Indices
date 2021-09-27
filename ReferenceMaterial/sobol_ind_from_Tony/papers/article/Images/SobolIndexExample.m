% Sobol Index Example

d=4;
%a=0.5*(1:d).^2; %the a_j decay more quickly
a=0.5*(1:d)/d; %the a_j decay less quickly

%% True answers for the test function
partialvariance = 1./(3.*(1+a).^2)
variance = prod(partialvariance+1)-1
sobol = partialvariance/variance

f = @(x) prod(bsxfun(@times,bsxfun(@plus,abs(4*x-2),a),1./(1+a)),2);
nval = zeros(1,d+2);
timeval = zeros(1,d+2);
[mu,out] = cubSobol_g(f,[zeros(1,d);ones(1,d)],'abstol',1e-4,'reltol',0);
mu
nval(1) = out.n;
timeval(1) = out.time;
[Efsq,out] = cubSobol_g(@(x) f(x).*f(x),[zeros(1,d);ones(1,d)], ...
   'abstol',1e-4,'reltol',0);
nval(2) = out.n;
timeval(2) = out.time;
varf = Efsq - mu^2
Emujsq = zeros(1,d);
for j = 1:d
   gj = @(x) f(x(:,1:d)).*f(x(:,[d+1:d+j-1 j d+j:2*d-1]));
   [Emujsq(j),out] = cubSobol_g(gj, ...
      [zeros(1,2*d-1);ones(1,2*d-1)],'abstol',1e-4,'reltol',0);
   nval(j+2) = out.n;
   timeval(j+2) = out.time;
end
varj = Efsq - Emujsq;
sobolind = (varf - varj)/varf
sobolerr = abs(sobol - sobolind)
nval
timeval
   

