function [q,int,out_param] = cubSobol_SI_g(varargin)
%CUBSOBOL_SI_G Quasi-Monte Carlo method using Sobol' cubatures
%to compute first order or total effect Sobol Indices
%within a specified generalized error tolerance with guarantees under
%Walsh-Fourier coefficients cone decay assumptions.
%
%   [q,out_param] = CUBSOBOL_SI_G(f,hyperbox,u) estimates the first order
%   Sobol Indices of f if u > 0, or the total effect Indice u otherwise,
%   where hyperbox is the sample space of the uniform distribution
%   (the distribution can also be set to normal), and the error guaranteed
%   not to be greater than a specific generalized error tolerance
%   tolfun:=max(abstol,reltol*| SI(f) |). Input f is a function handle. f should
%   accept an n x d matrix input, where d is the dimension and n is the 
%   number of points being evaluated simultaneously. The input hyperbox is
%   a 2 x d matrix, where the first row corresponds to the lower limits 
%   and the second row corresponds to the upper limits of the integral.
%   Given the construction of Sobol' sequences, d must be a positive 
%   integer with 1<=d<=1111.
%
%   q = CUBSOBOL_SI_G(f,hyperbox,u,measure,abstol,reltol)
%   estimates the u Sobol Indice of f. The answer is given within the
%   generalized error tolerance tolfun. All parameters should be
%   input in the order specified above. If an input is not specified,
%   the default value is used. Note that if an input is not specified,
%   the remaining tail cannot be specified either. Inputs f and hyperbox 
%   are required. The other optional inputs are in the correct order:
%   measure,abstol,reltol,mmin,mmax,fudge,toltype and
%   theta.
%
%   q = CUBSOBOL_SI_G(f,hyperbox,u,'measure',measure,'abstol',abstol,'reltol',reltol)
%   estimates the u Sobol Indice of f. The answer is given within
%   the generalized error tolerance tolfun. All the field-value pairs
%   are optional and can be supplied in any order. If an input is not
%   specified, the default value is used.
%
%   q = CUBSOBOL_SI_G(f,hyperbox,u,in_param) estimates the u Sobol Indice
%   of f. The answer is given within the generalized error tolerance tolfun.
% 
%   Input Arguments
%
%     f --- the integrand whose input should be a matrix n x d where n is
%     the number of data points and d the dimension, which cannot be
%     greater than 1111. By default f is f=@ x.^2.
%
%     hyperbox --- sample space of the distribution that defines the random
%     input vector. It must be be a 2 x d matrix, where the first row
%     corresponds to the lower limits and the second row to the upper
%     limits. The default value is [0;1].
%
%     u --- Integer that defines the dimensions to be checked. It needs to
%     be different than 0 and the absolute value of u needs to be smaller
%     or equal to the dimension of the function. Negative values of u mean
%     total effect Indice (all dimensions except u). By default is 1.
%
%     in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%     measure of a uniformly distributed random variable in the hyperbox
%     (each dimension is independent) or normally distributed
%     with covariance matrix I_d. The only possible values are
%     'uniform' or 'normal'. For 'uniform', the hyperbox must have
%     a finite volume while for 'normal', the hyperbox can only be defined as 
%     (-Inf,Inf)^d. By default it is 'uniform'.
%
%     in_param.abstol --- the absolute error tolerance, abstol>=0. By 
%     default it is 1e-4.
%
%     in_param.reltol --- the relative error tolerance, which should be
%     in [0,1]. Default value is 1e-2.
% 
%   Optional Input Arguments
% 
%     in_param.mmin --- the minimum number of points to start is 2^mmin.
%     The cone condition on the Fourier coefficients decay requires a
%     minimum number of points to start. The advice is to consider at least
%     mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
%     default it is 10.
%
%     in_param.mmax --- the maximum budget is 2^mmax. By construction of
%     the Sobol' generator, mmax is a positive integer such that
%     mmin<=mmax<=53. The default value is 24.
%
%     in_param.fudge --- the positive function multiplying the finite 
%     sum of Fast Walsh Fourier coefficients specified in the cone of functions.
%     This input is a function handle. The fudge should accept an array of
%     nonnegative integers being evaluated simultaneously. For more
%     technical information about this parameter, refer to the references.
%     By default it is @(m) 5*2.^-m.
%
%     in_param.toltype --- this is the generalized tolerance function.
%     There are two choices, 'max' which takes
%     max(abstol,reltol*| integral(f) | ) and 'comb' which is the linear combination
%     theta*abstol+(1-theta)*reltol*| integral(f) | . Theta is another 
%     parameter to be specified with 'comb'(see below). For pure absolute
%     error, either choose 'max' and set reltol = 0 or choose 'comb' and set
%     theta = 1. For pure relative error, either choose 'max' and set 
%     abstol = 0 or choose 'comb' and set theta = 0. Note that with 'max',
%     the user can not input abstol = reltol = 0 and with 'comb', if theta = 1
%     abstol con not be 0 while if theta = 0, reltol can not be 0.
%     By default toltype is 'max'.
% 
%     in_param.theta --- this input is parametrizing the toltype 
%     'comb'. Thus, it is only active when the toltype
%     chosen is 'comb'. It establishes the linear combination weight
%     between the absolute and relative tolerances
%     theta*abstol+(1-theta)*reltol*| integral(f) |. Note that for theta = 1, 
%     we have pure absolute tolerance while for theta = 0, we have pure 
%     relative tolerance. By default, theta=1.
%
%   Output Arguments
%
%     q --- the estimated value of the Sobol Indice.
%
%     out_param.d --- dimension of the domain of f.
%
%     out_param.n --- number of Sobol' points used to compute the
%     Sobol Indice of f.
%
%     out_param.bound_err --- predicted error bound of q based on the cone
%     conditions. If the function lies in the cone, the real error will be
%     smaller than generalized tolerance.
%
%     out_param.time --- time elapsed in seconds when calling cubSobol_SI_g.
%
%     out_param.exitflag --- this is a binary vector stating whether
%     warning flags arise. These flags tell about which conditions make the
%     final result certainly not guaranteed. One flag is considered arisen
%     when its value is 1. The following list explains the flags in the
%     respective vector order:
%
%                       1 : If reaching overbudget. It states whether
%                       the max budget is attained without reaching the
%                       guaranteed error tolerance.
%      
%                       2 : If the function lies outside the cone. In
%                       this case, results are not guaranteed. For more
%                       information about the cone definition, check the
%                       article mentioned below.
%
%     out_param.small --- Boolean indicating if we changed our estimator
%     for small first order Sobol Indices. This improves the estimation of
%     the Indices when they are small.
% 
%  Guarantee
% This algorithm computes first order and total effect Sobol Indices of
% real valued functions in [0,1)^d with a prescribed generalized error
% tolerance. The Walsh-Fourier coefficients of the integrand are assumed
% to be absolutely convergent. If the algorithm terminates without warning
% messages, the output is given with guarantees under the assumption that
% the integrand lies inside a cone of functions. The guarantee is based on 
% the decay rate of the Walsh-Fourier coefficients. For more details on how
% the cone is defined, please refer to the references below.
% 
%  Examples
% 
% Example 1:
% Ishigami example, first order indice for dimension 1:
% 
% >> f = @(x) sin(x(:,1)).*(1+1/10*(x(:,3).^4))+7*sin(x(:,2)).^2; hyperbox = pi*[-ones(1,3) ; ones(1,3)];
% >> q = cubSobol_SI_g(f,hyperbox,1,'uniform',1e-2,0); exactsol = .3139051827;
% >> check = abs(exactsol-q) < gail.tolfun(1e-2,0,1,exactsol,'max')
% check = 1
% 
% 
% Example 2:
% Bratley example, total effect indice for dimension 1:
% 
% >> f = @(x) sum(bsxfun(@times, cumprod(x,2), (-1).^(1:6)),2); hyperbox = [-zeros(1,6) ; ones(1,6)];
% >> q = cubSobol_SI_g(f,hyperbox,-1,'uniform',1e-2,1e-1); exactsol = .7396477462;
% >> check = abs(exactsol-q) < gail.tolfun(1e-2,1e-1,1,exactsol,'max')
% check = 1
% 
% 
% Example 3: 
% Sobol' g-function example, first order indice for dimension 5:
% 
% >> a = [0 1/2 3 9 99 99]; hyperbox = [zeros(1,6) ; ones(1,6)];
% >> f = @(x) prod(bsxfun(@(x,a) (abs(4*x-2)+a)./(1+a), x , a),2);
% >> q = cubSobol_SI_g(f,hyperbox,5,'uniform',1e-3,0); exactsol = 0.5867753221e-4;
% >> check = abs(exactsol-q) < gail.tolfun(1e-3,0,1,exactsol,'max')
% check = 1
%
%
% Example 4: 
% Morokoff and Caflish example, total effect indice for dimension 4:
% 
% >> f = @(x) (1+1/6)^6*prod(x,2).^(1/6); hyperbox = [zeros(1,6) ; ones(1,6)];
% >> q = cubSobol_SI_g(f,hyperbox,-4,'uniform',1e-3,1e-2); exactsol = .1753745708;
% >> check = abs(exactsol-q) < gail.tolfun(1e-3,1e-2,1,exactsol,'max')
% check = 1
%
%
%   See also CUBSOBOL_G
% 
%  References
%
%   [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama, "Reliable adaptive
%   cubature using digital sequences," 2014. Submitted for publication:
%   arXiv:1410.8615.
%
%   [2] Art B. Owen, "Better Estimation of Small Sobol' Sensitivity
%   Indices," ACM Trans. Model. Comput. Simul., 23, 2, Article 11 (May 2013).
%
%   [3] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   GAIL: Guaranteed Automatic Integration Library (Version 2.1)
%   [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
%   [4] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [5] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://code.google.com/p/gail/ 
%
%   [6] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
%   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
%   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
%   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
%   Workshop On Sustainable Software for Science: Practice And Experiences
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%

t_start = tic;
%% Initial important cone factors and Check-initialize parameters
r_lag = 4; %distance between coefficients summed and those computed
threshold_small = 0.1; % Threshold for changing the estimators for small indices
[f,hyperbox,out_param] = cubSobol_SI_g_param(r_lag,varargin{:});
u = out_param.u;
converged = false(2, out_param.d); % We flag the indices that converged
l_star = out_param.mmin - r_lag; % Minimum gathering of points for the sums of DFWT
omg_circ = @(m) 2.^(-m);
omg_hat = @(m) out_param.fudge(m)/((1+out_param.fudge(r_lag))*omg_circ(r_lag));

if strcmp(out_param.measure,'normal')
   f=@(x) f(gail.stdnorminv(x));
elseif strcmp(out_param.measure,'uniform')
   Cnorm = prod(hyperbox(2,:)-hyperbox(1,:));
   f=@(x) Cnorm*f(bsxfun(@plus,hyperbox(1,:),bsxfun(@times,(hyperbox(2,:)-hyperbox(1,:)),x))); % a + (b-a)x = u
end

% First and total order indices estimators and function evaluations
Sfo  = @(b,c) min(max(max(b(:,1),0)./max((c(:,2)-b(:,3).^2),0),0),1); % S function for first order
ffo = @(xpts,u,fx,fy,fxy) (fx - 1/2*mean(fx + fxy)).*(fxy - 1/2*mean(fx + fxy));
Sfo_s = Sfo; % S function for first order small
ffo_s = @(xpts,u,fx,fy,fxy) (fx - f([xpts(:,1:u-1) xpts(:,2*out_param.d + u) xpts(:,u+1:out_param.d)])).*(fxy - fy); % We redefine the non normalized estimator
Stot = Sfo; % S function for total effect
ftot = @(xpts,u,fx,fy,fxy) 1/2*(fy - fxy).^2;


%% Main algorithm (we added 2xd dimensions for each index and and additional 2 for mean of f and f^2)
sobstr = sobolset(3*out_param.d); %generate a Sobol' sequence 3*d to consider the changing to the estimator for some smaller size indices
sobstr = scramble(sobstr,'MatousekAffineOwen'); %scramble it
kappanumap_fx2_fx = bsxfun(@times,(1:2^out_param.mmin)', [1 1]); %initialize map
Stilde_fx2_fx = zeros(out_param.mmax-out_param.mmin+1, 2); %initialize sum of DFWT terms for fx2 and fx
CStilde_low_fx2_fx = -inf(out_param.mmax-l_star+1, 2); %initialize various sums of DFWT terms for necessary conditions for fx2 and fx
CStilde_up_fx2_fx = inf(out_param.mmax-l_star+1, 2); %initialize various sums of DFWT terms for necessary conditions for fx2 and fx
err_bound_int_fx2_fx = inf(out_param.mmax-out_param.mmin+1, 2); %initialize error estimates for fx2 and fx
est_int_fx2_fx = zeros(out_param.mmax-out_param.mmin+1, 2); % Estimate of the mean for fx2 and fx
exit_len = 2;
out_param.exit = false(exit_len,1); %we start the algorithm with all warning flags down
out_param.small = false; % boolean indicating if use the small indice estimator or not

INDICE.kappanumap = kappanumap_fx2_fx(:,1);
INDICE.Stilde = zeros(out_param.mmax-out_param.mmin+1,1);
INDICE.CStilde_low = CStilde_low_fx2_fx(:,1);
INDICE.CStilde_up = CStilde_up_fx2_fx(:,1);
if u > 0
    INDICE.S = Sfo;
    INDICE.f = ffo;
else
    INDICE.S = Stot;
    INDICE.f = ftot;
end
INDICE.est_int = est_int_fx2_fx(:,1); %initialize mean estimates for the integral in the numerator
INDICE.err_bound_int = err_bound_int_fx2_fx(:,1); %initialize error estimates for the integral in the numerator
INDICE.errest = inf(out_param.mmax-out_param.mmin+1); %initialize error estimates for the indice


%% Initial points and FWT
out_param.n = 2^out_param.mmin; %total number of points to start with
n0 = out_param.n; %initial number of points
xpts = sobstr(1:n0,1:3*out_param.d); %grab Sobol' points
fx = f(xpts(:,1:out_param.d)); %evaluate integrands y3
fx2 = fx.^2; %evaluate integrands y2
fxval = fx; % We store fx because fx will later contain the fast transform
fx2val = fx2; % We store fx2 because fx2 will later contain the fast transform
fy = f(xpts(:,out_param.d+1:2*out_param.d)); % We evaluate the f at the replicated points


fxy = f([xpts(:,out_param.d+1:out_param.d+abs(u)-1) xpts(:,abs(u)) xpts(:,out_param.d+abs(u)+1:2*out_param.d)]);
INDICE.y = INDICE.f(xpts,u,fx,fy,fxy);
if u > 0 % Check if the index is small
    aux_double = INDICE.S([mean(INDICE.y) mean(fx2) mean(fx)],[mean(INDICE.y) mean(fx2) mean(fx)]);
    if aux_double < threshold_small % If the normalized first order index is small, we change to a better estimator estimator
        INDICE.S = Sfo_s; % We redefine the S function for the small estimator
        INDICE.f = ffo_s;
        INDICE.y = INDICE.f(xpts,u,fx,fy,fxy); % We reevaluate the points if we change the estimator
        out_param.small = true;
    end
end
INDICE.est_int(1) = mean(INDICE.y); % Estimate the integral


%% Compute initial FWT
nllstart = int64(2^(out_param.mmin-r_lag-1));
for l=0:out_param.mmin-1 % We need the FWT for fx, fx^2, and the y values (a1/[a2-a3])
    nl=2^l;
    nmminlm1=2^(out_param.mmin-l-1);
    ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
    evenval=INDICE.y(ptind);
    oddval=INDICE.y(~ptind);
    INDICE.y(ptind)=(evenval+oddval)/2;
    INDICE.y(~ptind)=(evenval-oddval)/2;
    evenval=fx(ptind);
    oddval=fx(~ptind);
    fx(ptind)=(evenval+oddval)/2;
    fx(~ptind)=(evenval-oddval)/2;
    evenval=fx2(ptind);
    oddval=fx2(~ptind);
    fx2(ptind)=(evenval+oddval)/2;
    fx2(~ptind)=(evenval-oddval)/2;
end
%y now contains the FWT coefficients 

%% Create kappanumap implicitly from the data
for l=out_param.mmin-1:-1:1
   nl=2^l;
   oldone=abs(INDICE.y(INDICE.kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(INDICE.y(INDICE.kappanumap(nl+2:2*nl))); %later values of kappa
   flip=find(newone>oldone); %which in the pair are the larger ones
   if ~isempty(flip)
       flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
       temp = INDICE.kappanumap(nl+1+flipall); %then flip 
       INDICE.kappanumap(nl+1+flipall) = INDICE.kappanumap(1+flipall); %them
       INDICE.kappanumap(1+flipall) = temp; %around
   end
   oldone=abs(fx2(kappanumap_fx2_fx(2:nl,1))); %earlier values of kappa, don't touch first one
   newone=abs(fx2(kappanumap_fx2_fx(nl+2:2*nl,1))); %later values of kappa
   flip=find(newone>oldone); %which in the pair are the larger ones
   if ~isempty(flip)
       flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
       temp = kappanumap_fx2_fx(nl+1+flipall,1); %then flip 
       kappanumap_fx2_fx(nl+1+flipall,1) = kappanumap_fx2_fx(1+flipall,1); %them
       kappanumap_fx2_fx(1+flipall,1) = temp; %around
   end
   oldone=abs(fx(kappanumap_fx2_fx(2:nl,2))); %earlier values of kappa, don't touch first one
   newone=abs(fx(kappanumap_fx2_fx(nl+2:2*nl,2))); %later values of kappa
   flip=find(newone>oldone); %which in the pair are the larger ones
   if ~isempty(flip)
       flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
       temp=kappanumap_fx2_fx(nl+1+flipall,2); %then flip 
       kappanumap_fx2_fx(nl+1+flipall,2) = kappanumap_fx2_fx(1+flipall,2); %them
       kappanumap_fx2_fx(1+flipall,2) = temp; %around
   end
end

%% Compute Stilde
% We keep the error estimates for int fx2
Stilde_fx2_fx(1,1) = sum(abs(fx2(kappanumap_fx2_fx(nllstart+1:2*nllstart,1))));
est_int_fx2_fx(1,1,end) = mean(fx2val); % Estimate the integral
err_bound_int_fx2_fx(1,1) = out_param.fudge(out_param.mmin)*Stilde_fx2_fx(1,1);
% We keep the error estimates for int fx
Stilde_fx2_fx(1,2) = sum(abs(fx(kappanumap_fx2_fx(nllstart+1:2*nllstart,2))));
est_int_fx2_fx(1,2) = mean(fxval); % Estimate the integral
err_bound_int_fx2_fx(1,2) = out_param.fudge(out_param.mmin)*Stilde_fx2_fx(1,2);
int = est_int_fx2_fx(1,2); % Estimate of the expectation of the function

INDICE.Stilde(1) = sum(abs(INDICE.y(INDICE.kappanumap(nllstart+1:2*nllstart))));
INDICE.err_bound_int(1) = out_param.fudge(out_param.mmin)*INDICE.Stilde(1);

b = [INDICE.est_int(1) - INDICE.err_bound_int(1), est_int_fx2_fx(1,1) - err_bound_int_fx2_fx(1,1), est_int_fx2_fx(1,2) - err_bound_int_fx2_fx(1,2)];
c = [INDICE.est_int(1) + INDICE.err_bound_int(1), est_int_fx2_fx(1,1) + err_bound_int_fx2_fx(1,1), est_int_fx2_fx(1,2) + err_bound_int_fx2_fx(1,2)];
q = 1/2*(min(INDICE.S(c,b),1) + max(INDICE.S(b,c),0));
out_param.bound_err = 1/2*(min(INDICE.S(c,b),1) - max(INDICE.S(b,c),0));
INDICE.errest(1) = out_param.bound_err;


% Necessary conditions for all indices integrals and fx and fx2
for l = l_star:out_param.mmin % Storing the information for the necessary conditions
    C_low = 1/(1+omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    C_up = 1/(1-omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    INDICE.CStilde_low(l-l_star+1) = max(INDICE.CStilde_low(l-l_star+1),C_low*sum(abs(INDICE.y(INDICE.kappanumap(2^(l-1)+1:2^l)))));
    if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
        INDICE.CStilde_up(l-l_star+1) = min(INDICE.CStilde_up(l-l_star+1),C_up*sum(abs(INDICE.y(INDICE.kappanumap(2^(l-1)+1:2^l)))));
    end
    CStilde_low_fx2_fx(l-l_star+1,1) = max(CStilde_low_fx2_fx(l-l_star+1,1),C_low*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
    CStilde_low_fx2_fx(l-l_star+1,2) = max(CStilde_low_fx2_fx(l-l_star+1,2),C_low*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
    if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
        CStilde_up_fx2_fx(l-l_star+1,1) = min(CStilde_up_fx2_fx(l-l_star+1,1),C_up*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
        CStilde_up_fx2_fx(l-l_star+1,2) = min(CStilde_up_fx2_fx(l-l_star+1,2),C_up*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
    end
end
aux_bool = any(CStilde_low_fx2_fx(:) > CStilde_up_fx2_fx(:)); % Variable that checks conditions violated for fx2 and fx
if any(INDICE.CStilde_low(:) > INDICE.CStilde_up(:)) || aux_bool
    out_param.exit(2) = true;
end

% Check the end of the algorithm
deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-INDICE.errest(1)),...
    out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+INDICE.errest(1)),out_param.toltype));
deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-INDICE.errest(1)),...
    out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+INDICE.errest(1)),out_param.toltype));

q = q+deltaminus;
if out_param.bound_err <= deltaplus
   converged = true;
elseif out_param.mmin == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
   out_param.exit(1) = true;
end
        
out_param.time=toc(t_start);

%% Loop over m
for m = out_param.mmin+1:out_param.mmax
    if converged
       break;
    end
    mnext=m-1;
    nnext=2^mnext;
    xnext=sobstr(n0+(1:nnext),1:3*out_param.d);
    n0=n0+nnext;
    meff=m-out_param.mmin+1;
    fxnext = f(xnext(:,1:out_param.d)); %evaluate integrands y3
    fy = f(xnext(:,out_param.d+1:2*out_param.d)); % We evaluate the f at the replicated points
    fx2next = fxnext.^2; %evaluate integrands y2
    fxy = f([xnext(:,out_param.d+1:out_param.d+abs(u)-1) xnext(:,abs(u)) xnext(:,out_param.d+abs(u)+1:2*out_param.d)]);
    out_param.n = 2^m;
    ynext = INDICE.f(xnext,u,fxnext,fy,fxy);
    INDICE.est_int(meff) = 1/2*(INDICE.est_int(meff-1) + mean(ynext)); % Estimate the integral
    est_int_fx2_fx(meff,1) = 1/2*(est_int_fx2_fx(meff-1,1) + mean(fx2next)); % Estimate the integral of f^2
    est_int_fx2_fx(meff,2) = 1/2*(est_int_fx2_fx(meff-1,2) + mean(fxnext)); % Estimate the integral of f
    int = est_int_fx2_fx(meff,2); % Estimate of the expectation of the function
    
    
    %% Compute initial FWT on next points
    nllstart=int64(2^(m-r_lag-1));
    for l=0:mnext-1
        nl=2^l;
        nmminlm1=2^(mnext-l-1);
        ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
        evenval=ynext(ptind);
        oddval=ynext(~ptind);
        ynext(ptind)=(evenval+oddval)/2;
        ynext(~ptind)=(evenval-oddval)/2;
        evenval=fxnext(ptind);
        oddval=fxnext(~ptind);
        fxnext(ptind)=(evenval+oddval)/2;
        fxnext(~ptind)=(evenval-oddval)/2;
        evenval=fx2next(ptind);
        oddval=fx2next(~ptind);
        fx2next(ptind)=(evenval+oddval)/2;
        fx2next(~ptind)=(evenval-oddval)/2;
    end
    %% Compute FWT on all points
    fx = [fx; fxnext];
    fx2 = [fx2; fx2next];
    INDICE.y = [INDICE.y; ynext];
    nl=2^mnext;
    ptind=[true(nl,1); false(nl,1)];
    evenval=INDICE.y(ptind);
    oddval=INDICE.y(~ptind);
    INDICE.y(ptind)=(evenval+oddval)/2;
    INDICE.y(~ptind)=(evenval-oddval)/2;
    evenval=fx(ptind);
    oddval=fx(~ptind);
    fx(ptind)=(evenval+oddval)/2;
    fx(~ptind)=(evenval-oddval)/2;
    evenval=fx2(ptind);
    oddval=fx2(~ptind);
    fx2(ptind)=(evenval+oddval)/2;
    fx2(~ptind)=(evenval-oddval)/2;
    
    %% Update kappanumap
    kappanumap_fx2_fx = [kappanumap_fx2_fx ; 2^(m-1) + kappanumap_fx2_fx]; %initialize map
    INDICE.kappanumap = [INDICE.kappanumap ; 2^(m-1) + INDICE.kappanumap]; %initialize map only for fx and fx2
    for l=m-1:-1:m-r_lag
        nl=2^l;
        oldone=abs(INDICE.y(INDICE.kappanumap(2:nl))); %earlier values of kappa, don't touch first one
        newone=abs(INDICE.y(INDICE.kappanumap(nl+2:2*nl))); %later values of kappa
        flip=find(newone>oldone);
        if ~isempty(flip)
            flipall=bsxfun(@plus,flip,0:2^(l+1):2^m-1);
            temp=INDICE.kappanumap(nl+1+flipall);
            INDICE.kappanumap(nl+1+flipall)=INDICE.kappanumap(1+flipall);
            INDICE.kappanumap(1+flipall)=temp;
        end
        oldone=abs(fx2(kappanumap_fx2_fx(2:nl,1))); %earlier values of kappa, don't touch first one
        newone=abs(fx2(kappanumap_fx2_fx(nl+2:2*nl,1))); %later values of kappa
        flip=find(newone>oldone);
        if ~isempty(flip)
            flipall=bsxfun(@plus,flip,0:2^(l+1):2^m-1);
            temp=kappanumap_fx2_fx(nl+1+flipall,1);
            kappanumap_fx2_fx(nl+1+flipall,1)=kappanumap_fx2_fx(1+flipall,1);
            kappanumap_fx2_fx(1+flipall,1)=temp;
        end
        oldone=abs(fx(kappanumap_fx2_fx(2:nl,2))); %earlier values of kappa, don't touch first one
        newone=abs(fx(kappanumap_fx2_fx(nl+2:2*nl,2))); %later values of kappa
        flip=find(newone>oldone);
        if ~isempty(flip)
            flipall=bsxfun(@plus,flip,0:2^(l+1):2^m-1);
            temp=kappanumap_fx2_fx(nl+1+flipall,2);
            kappanumap_fx2_fx(nl+1+flipall,2)=kappanumap_fx2_fx(1+flipall,2);
            kappanumap_fx2_fx(1+flipall,2)=temp;
        end
    end

    %% Compute Stilde
    % We keep the error estimates and integrals only for int fx2
    Stilde_fx2_fx(meff,1) = sum(abs(fx2(kappanumap_fx2_fx(nllstart+1:2*nllstart,1))));
    err_bound_int_fx2_fx(meff,1) = out_param.fudge(m)*Stilde_fx2_fx(meff,1);
    % We keep the error estimates and integrals only for int fx
    Stilde_fx2_fx(meff,2) = sum(abs(fx(kappanumap_fx2_fx(nllstart+1:2*nllstart,2))));
    err_bound_int_fx2_fx(meff,2) = out_param.fudge(m)*Stilde_fx2_fx(meff,2);
  
    INDICE.Stilde(meff) = sum(abs(INDICE.y(INDICE.kappanumap(nllstart+1:2*nllstart))));
    INDICE.err_bound_int(meff) = out_param.fudge(m)*INDICE.Stilde(meff); % Only error bound for the integral on the numerator

    b = [INDICE.est_int(meff) - INDICE.err_bound_int(meff), est_int_fx2_fx(meff,1) - err_bound_int_fx2_fx(meff,1), est_int_fx2_fx(meff,2) - err_bound_int_fx2_fx(meff,2)];
    c = [INDICE.est_int(meff) + INDICE.err_bound_int(meff), est_int_fx2_fx(meff,1) + err_bound_int_fx2_fx(meff,1), est_int_fx2_fx(meff,2) + err_bound_int_fx2_fx(meff,2)];

    q = 1/2*(min(INDICE.S(c,b),1) + max(INDICE.S(b,c),0));
    out_param.bound_err = 1/2*(min(INDICE.S(c,b),1) - max(INDICE.S(b,c),0));
    INDICE.errest(meff) = out_param.bound_err;


    % Necessary conditions for all indices integrals and fx and fx2
    for l = l_star:m % Storing the information for the necessary conditions
        C_low = 1/(1+omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
        C_up = 1/(1-omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
        INDICE.CStilde_low(l-l_star+1) = max(INDICE.CStilde_low(l-l_star+1),C_low*sum(abs(INDICE.y(INDICE.kappanumap(2^(l-1)+1:2^l)))));
        if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
            INDICE.CStilde_up(l-l_star+1) = min(INDICE.CStilde_up(l-l_star+1),C_up*sum(abs(INDICE.y(INDICE.kappanumap(2^(l-1)+1:2^l)))));
        end
        CStilde_low_fx2_fx(l-l_star+1,1) = max(CStilde_low_fx2_fx(l-l_star+1,1),C_low*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
        CStilde_low_fx2_fx(l-l_star+1,2) = max(CStilde_low_fx2_fx(l-l_star+1,2),C_low*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
        if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
            CStilde_up_fx2_fx(l-l_star+1,1) = min(CStilde_up_fx2_fx(l-l_star+1,1),C_up*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
            CStilde_up_fx2_fx(l-l_star+1,2) = min(CStilde_up_fx2_fx(l-l_star+1,2),C_up*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
        end
    end
    aux_bool = any(CStilde_low_fx2_fx(:) > CStilde_up_fx2_fx(:)); % Variable that checks conditions violated for fx2 and fx
    if (any(INDICE.CStilde_low(:) > INDICE.CStilde_up(:)) || aux_bool)
        out_param.exit(2) = true;
    end

    % Check the end of the algorithm
    deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-INDICE.errest(meff)),...
        out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+INDICE.errest(meff)),out_param.toltype));
    deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-INDICE.errest(meff)),...
        out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+INDICE.errest(meff)),out_param.toltype));

    q = q + deltaminus;
    if out_param.bound_err <= deltaplus
       converged = true;
    elseif m == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
       out_param.exit(1) = true;
    end

    out_param.time=toc(t_start);
end

% Decode the exit structure
out_param.exitflag = (2.^(0:exit_len-1))'.*out_param.exit(:);
out_param = rmfield(out_param,'exit');

out_param.time=toc(t_start);
end


%% Parsing for the input of cubSobol_SI_all_g
function [f,hyperbox, out_param] = cubSobol_SI_g_param(r_lag,varargin)

% Default parameter values
default.u = 1;
default.hyperbox = [zeros(1,1);ones(1,1)];% default hyperbox
default.measure  = 'uniform';
default.abstol  = 1e-4;
default.reltol  = 1e-2;
default.mmin  = 10;
default.mmax  = 24;
default.fudge = @(m) 5*2.^-m;
default.toltype  = 'max';
default.theta  = 1;

if numel(varargin)<3
    help cubSobol_SI_g
    warning('GAIL:cubSobol_SI_g:fdnotgiven',...
        'At least, dimension u, function f and hyperbox need to be specified. Example for f(x) = x^2 and u = 1:')
    f = @(x) x.^2;
    out_param.f=f;
    out_param.u = default.u;
    hyperbox = default.hyperbox;
else
    f = varargin{1};
    if ~gail.isfcn(f)
        warning('GAIL:cubSobol_SI_g:fnotfcn',...
            'The given input f was not a function. Example for f(x) = x^2 and u = 1:')
        f = @(x) x.^2;
        out_param.f=f;
        out_param.u = default.u;
        hyperbox = default.hyperbox;
    else
        out_param.f=f;
        hyperbox = varargin{2};
        if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(size(hyperbox,2)<1111)
            warning('GAIL:cubSobol_SI_g:hyperbox_error1',...
                'The hyperbox must be a real matrix of size 2xd where d can not be greater than 1111. Example for f(x) = x^2 and u = 1:')
            f = @(x) x.^2;
            out_param.f=f;
            out_param.u = default.u;
            hyperbox = default.hyperbox;
        else
            u = varargin{3};
            if ~(ceil(u)==u) || abs(u) > size(hyperbox,2) || u == 0
                warning('GAIL:cubSobol_SI_g:u_error1',...
                    'Dimension u must be a an integer such that -d <= u <= d and u != 0. Example for f(x) = x^2 and u = 1:')
                f = @(x) x.^2;
                out_param.f=f;
                out_param.u = default.u;
                hyperbox = default.hyperbox;
            end
        end
    end
end

validvarargin=numel(varargin)>3;
if validvarargin
    in4=varargin(4:end);
    for j=1:numel(varargin)-3
    validvarargin=validvarargin && (isnumeric(in4{j}) ...
        || ischar(in4{j}) || isstruct(in4{j}) || gail.isfcn(in4{j}));
    end
    if ~validvarargin
        warning('GAIL:cubSobol_SI_g:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in4=varargin{4};
end

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
  f_addParamVal = @addParameter;
else
  f_addParamVal = @addParamValue;
end

if ~validvarargin
    out_param.measure = default.measure;
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;  
    out_param.fudge = default.fudge;
    out_param.toltype = default.toltype;
    out_param.theta = default.theta;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'hyperbox',@isnumeric);
    addRequired(p,'u',@isnumeric);
    if isnumeric(in4) || ischar(in4)
        addOptional(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@gail.isfcn);
        addOptional(p,'toltype',default.toltype,...
            @(x) any(validatestring(x, {'max','comb'})));
        addOptional(p,'theta',default.theta,@isnumeric);
    else
        if isstruct(in4) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'reltol',default.reltol,@isnumeric);
        f_addParamVal(p,'mmin',default.mmin,@isnumeric);
        f_addParamVal(p,'mmax',default.mmax,@isnumeric);
        f_addParamVal(p,'fudge',default.fudge,@gail.isfcn);
        f_addParamVal(p,'toltype',default.toltype,...
            @(x) any(validatestring(x, {'max','comb'})));
        f_addParamVal(p,'theta',default.theta,@isnumeric);
    end
    parse(p,f,hyperbox,u,varargin{4:end})
    out_param = p.Results;
end

out_param.d = size(hyperbox,2);

fdgyes = 0; % We store how many functions are in varargin. There can only
            % two functions as input, the function f and the fudge factor.
for j = 1:size(varargin,2)
    fdgyes = gail.isfcn(varargin{j})+fdgyes;
end
if fdgyes < 2 % No fudge factor given as input
    default.fudge = @(m) 5*2.^-(m/d);
end

%hyperbox should be 2 x dimension
if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(out_param.d<1111)
    warning('GAIL:cubSobol_SI_g:hyperbox_error2',...
        'The hyperbox must be a real matrix of size 2 x d where d can not be greater than 1111. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    hyperbox = default.hyperbox;
end

% Force measure to be uniform or normal only
if ~(strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') )
    warning('GAIL:cubSobol_SI_g:notmeasure',['The measure can only be uniform or normal.' ...
            ' Using default measure ' num2str(default.measure)])
    out_param.measure = default.measure;
end

% Force absolute tolerance greater than 0
if (out_param.abstol < 0 )
    warning('GAIL:cubSobol_SI_g:abstolnonpos',['Absolute tolerance cannot be negative.' ...
            ' Using default absolute tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force relative tolerance greater than 0 and smaller than 1
if (out_param.reltol < 0) || (out_param.reltol > 1)
    warning('GAIL:cubSobol_SI_g:reltolnonunit',['Relative tolerance should be chosen in [0,1].' ...
            ' Using default relative tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('GAIL:cubSobol_SI_g:lowmmin',['The minimum starting exponent ' ...
            'should be an integer greater than 0 and smaller or equal than the maxium.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force mmin to be integer greater than r_lag (so that l_star=mmin-r_lag>=0)
if out_param.mmin < r_lag
    warning('GAIL:cubSobol_SI_g:lowmminrlag',['The minimum starting exponent ' ...
            'should be at least ' num2str(r_lag) '.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 54
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=53)
    warning('GAIL:cubSobol_SI_g:wrongmmax',['The maximum exponent for the budget should be an integer biger than mmin and smaller than 54.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('GAIL:cubSobol_SI_g:fudgenonpos',['The fudge factor should be a positive function.' ...
            ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end

% Force toltype to be max or comb
if ~(strcmp(out_param.toltype,'max') || strcmp(out_param.toltype,'comb') )
    warning('GAIL:cubSobol_SI_g:nottoltype',['The error type can only be max or comb.' ...
            ' Using default toltype ' num2str(default.toltype)])
    out_param.toltype = default.toltype;
end

% Force theta to be in [0,1]
if (out_param.theta < 0) || (out_param.theta > 1)
    warning('GAIL:cubSobol_SI_g:thetanonunit',['Theta should be chosen in [0,1].' ...
            ' Using default theta ' num2str(default.theta)])
    out_param.theta = default.theta;
end

% Checking on pure absolute/relative error
if (out_param.abstol==0) && (out_param.reltol==0)
    warning('GAIL:cubSobol_SI_g:tolzeros',['Absolute and relative error tolerances can not be simultaniusly 0.' ...
            ' Using default absolute tolerance ' num2str(default.abstol) ' and relative tolerance ' num2str(default.reltol)])
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
end
if (strcmp(out_param.toltype,'comb')) && (out_param.theta==1) && (out_param.abstol==0)
    warning('GAIL:cubSobol_SI_g:abstolzero',['When choosing toltype comb, if theta=1 then abstol>0.' ...
            ' Using default absolute tolerance ' num2str(default.abstol) ])
    out_param.abstol = default.abstol;
end
if (strcmp(out_param.toltype,'comb')) && (out_param.theta==0) && (out_param.reltol==0)
    warning('GAIL:cubSobol_SI_g:reltolzero',['When choosing toltype comb, if theta=0 then reltol>0.' ...
            ' Using default relative tolerance ' num2str(default.reltol) ])
    out_param.reltol = default.reltol;
end

% Checking on the hyperbox given the measure
if (strcmp(out_param.measure,'uniform')) && ~all(all(isfinite(hyperbox)))
    warning('GAIL:cubSobol_SI_g:hyperboxnotfinite',['If uniform measure, hyperbox must be of finite volume.' ...
            ' Using default hyperbox:'])
    disp([zeros(1,out_param.d);ones(1,out_param.d)])
    hyperbox = [zeros(1,out_param.d);ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (sum(sum(isfinite(hyperbox)))>0)
    warning('GAIL:cubSobol_SI_g:hyperboxfinite',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (any(hyperbox(1,:)==hyperbox(2,:)) || any(hyperbox(1,:)>hyperbox(2,:)))
    warning('GAIL:cubSobol_SI_g:hyperboxnormalwrong',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end

end
