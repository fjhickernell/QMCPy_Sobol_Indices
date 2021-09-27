function [q,int,out_param] = cubSobol_SI_fo_g(varargin)
%CUBSOBOL_SI_FO_G Quasi-Monte Carlo method using Sobol' cubatures
%to compute all first order Sobol Indices whitin a specified
%generalized error tolerance with guarantees under Walsh-Fourier
%coefficients cone decay assumptions.
%
%   [q,out_param] = CUBSOBOL_SI_FO_G(f,hyperbox) estimates all first order
%   Sobol Indices of f, where hyperbox is the sample space
%   of the uniform distribution (the distribution can also be set to normal),
%   and the estimation error is guaranteed not to be greater than a specific
%   generalized error tolerance tolfun:=max(abstol,reltol*| SI(f) |).
%   Input f is a function handle. f should accept an n x d matrix input,
%   where d is the dimension and n is the number of points being evaluated
%   simultaneously. The input hyperbox is a 2 x d matrix, where the first
%   row corresponds to the lower limits and the second row corresponds to
%   the upper limits of the integral. Given the construction of Sobol'
%   sequences, d must be a positive integer with 1<=d<=370.
%
%   q = CUBSOBOL_SI_FO_G(f,hyperbox,measure,abstol,reltol) estimates all
%   first order Sobol Indices of f. The answer is given within
%   the generalized error tolerance tolfun. All parameters should be
%   input in the order specified above. If an input is not specified,
%   the default value is used. Note that if an input is not specified,
%   the remaining tail cannot be specified either. Inputs f and hyperbox 
%   are required. The other optional inputs are in the correct order:
%   measure,abstol,reltol,mmin,mmax,fudge,toltype and
%   theta.
%
%   q = CUBSOBOL_SI_FO_G(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%   estimates all first order Sobol Indices of f. The answer is given
%   within the generalized error tolerance tolfun. All the field-value
%   pairs are optional and can be supplied in any order. If an input is not
%   specified, the default value is used.
%
%   q = CUBSOBOL_SI_FO_G(f,hyperbox,in_param) estimates all first order
%   Sobol Indices of f. The answer is given within the generlized error
%   tolerance tolfun.
% 
%   Input Arguments
%
%     f --- the integrand whose input should be a matrix n x d where n is
%     the number of data points and d the dimension, which cannot be
%     greater than 370. By default f is f=@ x.^2.
%
%     hyperbox --- sample space of the distribution that defines the random
%     input vector. It must be be a 2 x d matrix, where the first row
%     corresponds to the lower limits and the second row to the upper
%     limits. The default value is [0;1].
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
%     q --- the estimated value of all Sobol Indices in matrix form.
%     Each column is the dimension for which we estimate the indice.
%
%     out_param.d --- dimension of the domain of f.
%
%     out_param.n --- number of Sobol' points used to compute q.
%
%     out_param.bound_err --- predicted error bound of q based on the cone
%     conditions. If the function lies in the cone, the real error will be
%     smaller than generalized tolerance.
%
%     out_param.time --- time elapsed in seconds when calling cubSobol_SI_fo_g.
%
%     out_param.exitflag --- for each q, this is a binary vector stating
%     whether warning flags arise. This is a triple array element with
%     flags stored in the first component, and indice specified by the
%     other two. These flags tell about which conditions make the
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
%     the indices when they are small.
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
% Ishigami example:
% 
% >> f = @(x) sin(x(:,1)).*(1+1/10*(x(:,3).^4))+7*sin(x(:,2)).^2; hyperbox = pi*[-ones(1,3) ; ones(1,3)];
% >> q = cubSobol_SI_fo_g(f,hyperbox,'uniform',1e-1,0); exactsol = [.3139051827, .4424111333, 0.];
% >> check = sum(sum(abs(exactsol-q) < gail.tolfun(1e-1,0,1,exactsol,'max')))
% check = 3
% 
% 
% Example 2:
% Bratley example:
% 
% >> f = @(x) sum(bsxfun(@times, cumprod(x,2), (-1).^(1:6)),2); hyperbox = [-zeros(1,6) ; ones(1,6)];
% >> q = cubSobol_SI_fo_g(f,hyperbox,'uniform',1e-1,1e-1); exactsol = [.6528636616, .1791303924, 0.3701041165e-1, 0.1332374820e-1, 0.1480416466e-2, 0.1480416466e-2];
% >> check = sum(sum(abs(exactsol-q) < gail.tolfun(1e-1,1e-1,1,exactsol,'max')))
% check = 6
% 
% 
% Example 3: 
% Sobol' g-function example, first order indice for dimension 5:
% 
% >> a = [0 1/2 3 9 99 99]; hyperbox = [zeros(1,6) ; ones(1,6)];
% >> f = @(x) prod(bsxfun(@(x,a) (abs(4*x-2)+a)./(1+a), x , a),2);
% >> q = cubSobol_SI_fo_g(f,hyperbox,'uniform',1e-2,1e-1); exactsol = [.5867811893, .2607916397, 0.3667382378e-1, 0.5867811312e-2, 0.5867753221e-4, 0.5867753221e-4];
% >> check = sum(sum(abs(exactsol-q) < gail.tolfun(1e-2,1e-1,1,exactsol,'max')))
% check = 6
%
%
% Example 4:
% Morokoff and Caflish example, total effect indice for dimension 4:
% 
% >> f = @(x) (1+1/6)^6*prod(x,2).^(1/6); hyperbox = [zeros(1,6) ; ones(1,6)];
% >> q = cubSobol_SI_fo_g(f,hyperbox,'uniform',1e-2,1e-1); exactsol = [.1581948744, .1581948744, .1581948744, .1581948744, .1581948744, .1581948744];
% >> check = sum(sum(abs(exactsol-q) < gail.tolfun(1e-2,1e-1,1,exactsol,'max')))
% check = 6
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
[f,hyperbox,out_param] = cubSobol_SI_fo_g_param(r_lag,varargin{:});
converged = false(1,out_param.d); % We flag the indices that converged
l_star = out_param.mmin - r_lag; % Minimum gathering of points for the sums of DFWT
omg_circ = @(m) 2.^(-m);
omg_hat = @(m) out_param.fudge(m)/((1+out_param.fudge(r_lag))*omg_circ(r_lag));

if strcmp(out_param.measure,'normal')
   f=@(x) f(gail.stdnorminv(x));
elseif strcmp(out_param.measure,'uniform')
   Cnorm = prod(hyperbox(2,:)-hyperbox(1,:));
   f=@(x) Cnorm*f(bsxfun(@plus,hyperbox(1,:),bsxfun(@times,(hyperbox(2,:)-hyperbox(1,:)),x))); % a + (b-a)x = u
end

% Below, both numerator and denominator need to be positive. In addition,
% the final value needs to be between 0 and 1.
numerator_size = 1;
% Sfo_min = @(down, up) estimate_min(down, up); % S function for first order
% Sfo_max = @(down, up) estimate_max(down, up); % S function for first order
% ffo1 = @(fx,fy,fxy) fx.*fxy;
% ffo2 = @(fx,fy,fxy) 1/2*(fx + fxy);
Sfo_min = @(down, up) estimate_min([down(1) 0 down(2:end)], [up(1) 0 up(2:end)]); % S function for first order
Sfo_max = @(down, up) estimate_max([down(1) 0 down(2:end)], [up(1) 0 up(2:end)]); % S function for first order
ffo = @(fx,fy,fxy) fx.*(fxy - fy);

%% Main algorithm (we added 2xd dimensions for each index and and additional 2 for mean of f and f^2)
sobstr = sobolset(2*out_param.d); %generate a Sobol' sequence 2*d
sobstr_scr = sobolset(out_param.d); % This will be the first scrambled sequence...
    % The replicated sequence will have to use the same scrambling
sobstr_scr = scramble(sobstr_scr,'MatousekAffineOwen'); %scramble it
kappanumap_fx2_fx = bsxfun(@times,(1:2^out_param.mmin)', [1 1]); %initialize map
Stilde_fx2_fx = zeros(out_param.mmax-out_param.mmin+1, 2); %initialize sum of DFWT terms for fx2 and fx
CStilde_low_fx2_fx = -inf(out_param.mmax-l_star+1, 2); %initialize various sums of DFWT terms for necessary conditions for fx2 and fx
CStilde_up_fx2_fx = inf(out_param.mmax-l_star+1, 2); %initialize various sums of DFWT terms for necessary conditions for fx2 and fx
err_bound_int_fx2_fx = inf(out_param.mmax-out_param.mmin+1, 2); %initialize error estimates for fx2 and fx
est_int_fx2_fx = zeros(out_param.mmax-out_param.mmin+1, 2); % Estimate of the mean for fx2 and fx
exit_len = 2;
out_param.exit = false(exit_len, out_param.d); %we start the algorithm with all warning flags down

for u = 1:out_param.d
    INDICES(u).kappanumap = kappanumap_fx2_fx(:, 1:numerator_size);
    INDICES(u).Stilde = zeros(out_param.mmax-out_param.mmin+1, numerator_size);
    INDICES(u).CStilde_low = CStilde_low_fx2_fx(:, 1:numerator_size);
    INDICES(u).CStilde_up = CStilde_up_fx2_fx(:, 1:numerator_size);
    INDICES(u).Smin = Sfo_min;
    INDICES(u).Smax = Sfo_max;
    if numerator_size == 1
        INDICES(u).f = {ffo}; % {ffo1, ffo2} {ffo}
    elseif numerator_size == 2
        INDICES(u).f = {ffo1, ffo2}; % {ffo1, ffo2} {ffo}
    else
        error('Check first order Sobol indices, number of functions in numerator')
    end
    INDICES(u).est_int = est_int_fx2_fx(:, 1:numerator_size); %initialize mean estimates for each integral in the numerator
    INDICES(u).err_bound_int = err_bound_int_fx2_fx(:, 1:numerator_size); %initialize error estimates for each integral in the numerator
    INDICES(u).errest = inf(out_param.mmax-out_param.mmin+1); %initialize error estimates for each indice
end


%% Initial points and FWT

%% In this section we apply the replicated method
% For the replication method, we want to sort dimension u+d with the same
% order as dimension u. In other words, if dimension u is sorted according
% to permutation tau and u+d according to sigma, we want to find psi such
% that sigma(psi) = tau. Therefore, psi = sigma^-1(tau). To find sigma^-1 we
% use the sorting Matlab function. The original order is sigma, so we only
% need to call the function values at order psi.
% Note that with scrambled points, to keep the replicated property,
% dimensions u and u+d need to share the same scrambling for all u (for
% each u we can choose a different scrambling). This means that the sorting
% psi is always the same that can be found with non-scr Sobol points (S for
% the Scrambling + shift): S(sigma(psi)) = S(tau) -> psi = 
% sigma^-1(S^-1(S(tau))) = sigma^-1(tau). Therefore, the sorting can be
% found with non-scr Sobol points and applied to the scrambled points.
out_param.n = 2^out_param.mmin*ones(1,out_param.d); %total number of points to start with
n0 = out_param.n(1); %initial number of points
xpts = sobstr(1:n0,1:2*out_param.d); %grab Sobol' points
tau = 2^(out_param.mmin)*xpts(:, 1:out_param.d) + 1;
sigma = 2^(out_param.mmin)*xpts(:, out_param.d+1:2*out_param.d) + 1;

%We create the replicated points
xpts_1 = sobstr_scr(1:n0, 1:out_param.d); xpts_2 = zeros(size(xpts_1));
for u = 1:out_param.d
    [~, tau_inv] = sort(tau(:, u));
    xpts_2(:, u) = xpts_1(tau_inv(sigma(:, u)), u);
end

% We perform some evaluations
fx = f(xpts_1); %evaluate integrands y3
fx2 = fx.^2; %evaluate integrands y2
fxval = fx; % We store fx because fx will later contain the fast transform
fx2val = fx2; % We store fx2 because fx2 will later contain the fast transform
% We only need fy if we use the estimator that includes ffo
fy = [];
if numerator_size == 1
    fy = f(xpts_2); % We evaluate the f at the replicated points
end

% We apply the replicated method
fxy_base = f(xpts_2); % In sigma order
for u = 1:out_param.d
    [~, sigma_inv] = sort(sigma(:, u));
    fxy = fxy_base(sigma_inv(tau(:, u)));
    INDICES(u).y = [];
    for p = 1:numerator_size
        INDICES(u).y = [INDICES(u).y INDICES(u).f{p}(fx,fy,fxy)];
    end
    INDICES(u).est_int(1, :) = mean(INDICES(u).y, 1); % Estimate the integral
end


%% Compute initial FWT
nllstart = int64(2^(out_param.mmin-r_lag-1));
for l=0:out_param.mmin-1 % We need the FWT for fx, fx^2, and the y values (a1/[a2-a3])
    nl=2^l;
    nmminlm1=2^(out_param.mmin-l-1);
    ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
    for u = 1:out_param.d
        evenval=INDICES(u).y(ptind, :);
        oddval=INDICES(u).y(~ptind, :);
        INDICES(u).y(ptind, :)=(evenval+oddval)/2;
        INDICES(u).y(~ptind, :)=(evenval-oddval)/2;
    end
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

for u = 1:out_param.d + 1 % u == d+1 is the kappanumap for fx2 and fx
    %% Create kappanumap implicitly from the data
    for l=out_param.mmin-1:-1:1
       nl=2^l;
       % for u = d+1 we compute the denominator which is the same for
       % all indices. That is why we store the y values and kappanumap
       % separetly for the denominator
       if u < out_param.d + 1
           for p = 1:numerator_size
                oldone=abs(INDICES(u).y(INDICES(u).kappanumap(2:nl, p), p)); %earlier values of kappa, don't touch first one
                newone=abs(INDICES(u).y(INDICES(u).kappanumap(nl+2:2*nl, p), p)); %later values of kappa
                flip=find(newone>oldone); %which in the pair are the larger ones
                if ~isempty(flip)
                   flipall=bsxfun(@plus,flip,0:2^(l+1):2^out_param.mmin-1);
                   temp = INDICES(u).kappanumap(nl+1+flipall, p); %then flip 
                   INDICES(u).kappanumap(nl+1+flipall, p) = INDICES(u).kappanumap(1+flipall, p); %them
                   INDICES(u).kappanumap(1+flipall, p) = temp; %around
                end
           end
       else % We keep the kappamap for int fx2 and fx
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

for u = 1:out_param.d
    INDICES(u).Stilde(1, :) = sum(abs(INDICES(u).y(INDICES(u).kappanumap(nllstart+1:2*nllstart, :))), 1);
    INDICES(u).err_bound_int(1, :) = out_param.fudge(out_param.mmin)*INDICES(u).Stilde(1, :);

    down = [INDICES(u).est_int(1, :) - INDICES(u).err_bound_int(1, :), est_int_fx2_fx(1, :) - err_bound_int_fx2_fx(1, :)];
    up = [INDICES(u).est_int(1, :) + INDICES(u).err_bound_int(1, :), est_int_fx2_fx(1, :) + err_bound_int_fx2_fx(1, :)];
    
    center = [INDICES(u).est_int(1, :), est_int_fx2_fx(1,1), est_int_fx2_fx(1,2)];
    out_param.center(u) = INDICES(u).Smax(center,center);
        
    q(u) = 1/2*(INDICES(u).Smax(down,up) + INDICES(u).Smin(down,up));
    out_param.bound_err(u) = 1/2*(INDICES(u).Smax(down,up) - INDICES(u).Smin(down,up));
    INDICES(u).errest(1) = out_param.bound_err(u);
end

% Necessary conditions for all indices integrals and fx and fx2
for l = l_star:out_param.mmin % Storing the information for the necessary conditions
    C_low = 1/(1+omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    C_up = 1/(1-omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
    for u = 1:out_param.d
        INDICES(u).CStilde_low(l-l_star+1, :) = max(INDICES(u).CStilde_low(l-l_star+1, :),C_low*sum(abs(INDICES(u).y(INDICES(u).kappanumap(2^(l-1)+1:2^l, :))), 1));
        if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
            INDICES(u).CStilde_up(l-l_star+1, :) = min(INDICES(u).CStilde_up(l-l_star+1, :),C_up*sum(abs(INDICES(u).y(INDICES(u).kappanumap(2^(l-1)+1:2^l, :))), 1));
        end
    end
    CStilde_low_fx2_fx(l-l_star+1,1) = max(CStilde_low_fx2_fx(l-l_star+1,1),C_low*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
    CStilde_low_fx2_fx(l-l_star+1,2) = max(CStilde_low_fx2_fx(l-l_star+1,2),C_low*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
    if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
        CStilde_up_fx2_fx(l-l_star+1,1) = min(CStilde_up_fx2_fx(l-l_star+1,1),C_up*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
        CStilde_up_fx2_fx(l-l_star+1,2) = min(CStilde_up_fx2_fx(l-l_star+1,2),C_up*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
    end
end
aux_bool = any(any(CStilde_low_fx2_fx > CStilde_up_fx2_fx)); % Variable that checks conditions violated for fx2 and fx
for u = 1:out_param.d
    if any(any(INDICES(u).CStilde_low > INDICES(u).CStilde_up)) || aux_bool
        out_param.exit(2,u) = true;
    end
end

for u = 1:out_param.d
    % Check the end of the algorithm
    deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q(u)-INDICES(u).errest(1)),...
        out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q(u)+INDICES(u).errest(1)),out_param.toltype));
    deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q(u)-INDICES(u).errest(1)),...
        out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q(u)+INDICES(u).errest(1)),out_param.toltype));

    q(u) = q(u)+deltaminus;
    if out_param.bound_err(u) <= deltaplus
       converged(u) = true;
    elseif out_param.mmin == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
       out_param.exit(1,u) = true;
    end
end
out_param.time=toc(t_start);

%% Loop over m
for m = out_param.mmin+1:out_param.mmax
    if all(converged)
       break;
    end
    mnext=m-1;
    nnext=2^mnext;
    xnext=sobstr(n0+(1:nnext),1:2*out_param.d); % The following Sobol 
                    % points have the nice property that can give us the
                    % order doing this operation. The +1 is because in
                    % Matlab, vectors start at 1, not 0.
    tau = (2^m*xnext(:, 1:out_param.d) - 1)/2 + 1;
    sigma = (2^m*xnext(:, out_param.d+1:2*out_param.d) - 1)/2 + 1;
    %We create the replicated points
    xnext_1 = sobstr_scr(n0+(1:nnext),1:out_param.d); xnext_2 = zeros(size(xnext_1));
    for u = 1:out_param.d
        [~, tau_inv] = sort(tau(:, u));
        xnext_2(:, u) = xnext_1(tau_inv(sigma(:, u)), u);
    end

    n0=n0+nnext;
    fxnext = f(xnext_1); %evaluate integrands y3
    if numerator_size == 1
        fy = f(xnext_2); % We evaluate the f at the replicated points
    end
    fx2next = fxnext.^2; %evaluate integrands y2
    fxval = [fxval; fxnext];
	fx2val = [fx2val; fx2next];
    
    %% Compute initial FWT on next points only for fx and fx2
    % We need to split fx2 and fx from the indices because the number of
    % data points used will vary depending on the indice but we will always
    % still need fx and fx2. Keeping altogehter requires too many if's.
    nllstart=int64(2^(m-r_lag-1));
    meff=m-out_param.mmin+1;
    for l=0:mnext-1
        nl=2^l;
        nmminlm1=2^(mnext-l-1);
        ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
        evenval=fxnext(ptind);
        oddval=fxnext(~ptind);
        fxnext(ptind)=(evenval+oddval)/2;
        fxnext(~ptind)=(evenval-oddval)/2;
        evenval=fx2next(ptind);
        oddval=fx2next(~ptind);
        fx2next(ptind)=(evenval+oddval)/2;
        fx2next(~ptind)=(evenval-oddval)/2;
    end
    %% Compute FWT on all points only for fx and fx2
    fx = [fx; fxnext];
    fx2 = [fx2; fx2next];
    nl=2^mnext;
    ptind=[true(nl,1); false(nl,1)];
    evenval=fx(ptind);
    oddval=fx(~ptind);
    fx(ptind)=(evenval+oddval)/2;
    fx(~ptind)=(evenval-oddval)/2;
    evenval=fx2(ptind);
    oddval=fx2(~ptind);
    fx2(ptind)=(evenval+oddval)/2;
    fx2(~ptind)=(evenval-oddval)/2;
    %% Update kappanumap only for fx and fx2
    kappanumap_fx2_fx = [kappanumap_fx2_fx ; 2^(m-1) + kappanumap_fx2_fx]; %initialize map only for fx and fx2
    for l=m-1:-1:m-r_lag
        nl=2^l;
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

    %% Compute Stilde for fx and fx2 only
    % We keep the error estimates and integrals only for int fx2
    Stilde_fx2_fx(meff,1) = sum(abs(fx2(kappanumap_fx2_fx(nllstart+1:2*nllstart,1))));
    est_int_fx2_fx(meff,1) = mean(fx2val); % Estimate the integral of f^2
    err_bound_int_fx2_fx(meff,1) = out_param.fudge(m)*Stilde_fx2_fx(meff,1);
    % We keep the error estimates and integrals only for int fx
    Stilde_fx2_fx(meff,2) = sum(abs(fx(kappanumap_fx2_fx(nllstart+1:2*nllstart,2))));
    est_int_fx2_fx(meff,2) = mean(fxval); % Estimate the integral of f
    err_bound_int_fx2_fx(meff,2) = out_param.fudge(m)*Stilde_fx2_fx(meff,2);
    int = est_int_fx2_fx(meff,2); % Estimate of the expectation of the function

    %% We start computing everything for the numerator. fx and fx2 were for
    % the denominator only, which is common with all indices
    fxy_base = f(xnext_2); % In sigma order,
        % and at this point, we now that at least 1 indice did not meet the
        % error tolerance so we need to evaluate fxy_base. However, the
        % sorting, FFWT, etc, only needs to be done for those that did not
        % converge.
    for u = 1:out_param.d
            if ~converged(u)
                INDICES(u).n = 2^m;
                out_param.n(u) = 2^m;
                [~, sigma_inv] = sort(sigma(:, u));
                fxy = fxy_base(sigma_inv(tau(:, u)));
                ynext = [];
                % fxval already was computed for fx and fx2 above, so fx
                % contains FFWT coefficients of fx.
                for p = 1:numerator_size
                    ynext = [ynext INDICES(u).f{p}(fxval(end/2 + 1:end),fy,fxy)];
                end
                INDICES(u).est_int(meff, :) = 1/2*(INDICES(u).est_int(meff-1, :) + mean(ynext, 1)); % Estimate the integral

               %% Compute initial FWT on next points only for y
               nllstart=int64(2^(m-r_lag-1));
               meff=m-out_param.mmin+1;
                for l=0:mnext-1
                    nl=2^l;
                    nmminlm1=2^(mnext-l-1);
                    ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
                    evenval=ynext(ptind, :);
                    oddval=ynext(~ptind, :);
                    ynext(ptind, :)=(evenval+oddval)/2;
                    ynext(~ptind, :)=(evenval-oddval)/2;
                end
                %% Compute FWT on all points only for y
                INDICES(u).y = [INDICES(u).y; ynext];
                nl=2^mnext;
                ptind=[true(nl,1); false(nl,1)];
                evenval = INDICES(u).y(ptind, :);
                oddval = INDICES(u).y(~ptind, :);
                INDICES(u).y(ptind, :) = (evenval+oddval)/2;
                INDICES(u).y(~ptind, :) = (evenval-oddval)/2;

                %% Update kappanumap only for indices
                INDICES(u).kappanumap = [INDICES(u).kappanumap ; 2^(m-1)+INDICES(u).kappanumap];
                for l=m-1:-1:m-r_lag
                  nl=2^l;
                  for p = 1:numerator_size
                      oldone=abs(INDICES(u).y(INDICES(u).kappanumap(2:nl, p), p)); %earlier values of kappa, don't touch first one
                      newone=abs(INDICES(u).y(INDICES(u).kappanumap(nl+2:2*nl, p), p)); %later values of kappa
                      flip = find(newone>oldone);
                      if ~isempty(flip)
                          flipall = bsxfun(@plus,flip,0:2^(l+1):2^m-1);
                          temp = INDICES(u).kappanumap(nl+1+flipall, p);
                          INDICES(u).kappanumap(nl+1+flipall, p) = INDICES(u).kappanumap(1+flipall, p);
                          INDICES(u).kappanumap(1+flipall, p) = temp;
                      end
                  end
                end
            
                %% Compute Stilde for y only
                INDICES(u).Stilde(meff, :) = sum(abs(INDICES(u).y(INDICES(u).kappanumap(nllstart+1:2*nllstart, :))), 1);
                INDICES(u).err_bound_int(meff, :) = out_param.fudge(m)*INDICES(u).Stilde(meff, :); % Only error bound for the integral on the numerator

                down = [INDICES(u).est_int(meff, :) - INDICES(u).err_bound_int(meff, :), est_int_fx2_fx(meff, :) - err_bound_int_fx2_fx(meff, :)];
                up = [INDICES(u).est_int(meff, :) + INDICES(u).err_bound_int(meff, :), est_int_fx2_fx(meff, :) + err_bound_int_fx2_fx(meff, :)];
                
                center = [INDICES(u).est_int(1, :), est_int_fx2_fx(1,1), est_int_fx2_fx(1,2)];
                out_param.center(u) = INDICES(u).Smax(center,center);
                
                q(u) = 1/2*(INDICES(u).Smax(down,up) + INDICES(u).Smin(down,up));
                out_param.bound_err(u) = 1/2*(INDICES(u).Smax(down,up) - INDICES(u).Smin(down,up));
                INDICES(u).errest(meff) = out_param.bound_err(u);
            end
    end

    % Necessary conditions for all indices integrals and fx and fx2
    for l = l_star:m % Storing the information for the necessary conditions
        C_low = 1/(1+omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
        C_up = 1/(1-omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l));
        for u = 1:out_param.d
            if ~converged(u)
                INDICES(u).CStilde_low(l-l_star+1, :) = max(INDICES(u).CStilde_low(l-l_star+1, :),C_low*sum(abs(INDICES(u).y(INDICES(u).kappanumap(2^(l-1)+1:2^l, :))), 1));
                if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
                    INDICES(u).CStilde_up(l-l_star+1, :) = min(INDICES(u).CStilde_up(l-l_star+1, :),C_up*sum(abs(INDICES(u).y(INDICES(u).kappanumap(2^(l-1)+1:2^l, :))), 1));
                end
            end
        end
        CStilde_low_fx2_fx(l-l_star+1,1) = max(CStilde_low_fx2_fx(l-l_star+1,1),C_low*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
        CStilde_low_fx2_fx(l-l_star+1,2) = max(CStilde_low_fx2_fx(l-l_star+1,2),C_low*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
        if (omg_hat(out_param.mmin-l)*omg_circ(out_param.mmin-l) < 1)
            CStilde_up_fx2_fx(l-l_star+1,1) = min(CStilde_up_fx2_fx(l-l_star+1,1),C_up*sum(abs(fx2(kappanumap_fx2_fx(2^(l-1)+1:2^l,1)))));
            CStilde_up_fx2_fx(l-l_star+1,2) = min(CStilde_up_fx2_fx(l-l_star+1,2),C_up*sum(abs(fx(kappanumap_fx2_fx(2^(l-1)+1:2^l,2)))));
        end
    end
    aux_bool = any(any(CStilde_low_fx2_fx > CStilde_up_fx2_fx)); % Variable that checks conditions violated for fx2 and fx
    for u = 1:out_param.d
        if ~converged(u) & (any(any(INDICES(u).CStilde_low > INDICES(u).CStilde_up)) || aux_bool)
            out_param.exit(2,u) = true;
        end
    end

    for u = 1:out_param.d
            if ~converged(u)
                % Check the end of the algorithm
                deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
                    out_param.reltol,out_param.theta,abs(q(u)-INDICES(u).errest(meff)),...
                    out_param.toltype)+gail.tolfun(out_param.abstol,out_param.reltol,...
                    out_param.theta,abs(q(u)+INDICES(u).errest(meff)),out_param.toltype));
                deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
                    out_param.reltol,out_param.theta,abs(q(u)-INDICES(u).errest(meff)),...
                    out_param.toltype)-gail.tolfun(out_param.abstol,out_param.reltol,...
                    out_param.theta,abs(q(u)+INDICES(u).errest(meff)),out_param.toltype));        

                q(u) = q(u) + deltaminus;
                if out_param.bound_err(u) <= deltaplus
                   converged(u) = true;
                elseif m == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
                   out_param.exit(1,u) = true;
                end
            end
    end
    
    out_param.time=toc(t_start);
end

% Decode the exit structure
for u = 1:out_param.d
        out_param.exitflag(:,u) = (2.^(0:exit_len-1))'.*out_param.exit(:,u);
end
out_param = rmfield(out_param,'exit');

out_param.time=toc(t_start);
end


%% Parsing for the input of cubSobol_SI_fo_g
function [f,hyperbox, out_param] = cubSobol_SI_fo_g_param(r_lag,varargin)

% Default parameter values
default.hyperbox = [zeros(1,1);ones(1,1)];% default hyperbox
default.measure  = 'uniform';
default.abstol  = 1e-4;
default.reltol  = 1e-2;
default.mmin  = 10;
default.mmax  = 24;
default.fudge = @(m) 5*2.^-m;
default.toltype  = 'max';
default.theta  = 1;

if numel(varargin)<2
    help cubSobol_SI_fo_g
    warning('GAIL:cubSobol_SI_fo_g:fdnotgiven',...
        'At least, function f and hyperbox need to be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    hyperbox = default.hyperbox;
else
    f = varargin{1};
    if ~gail.isfcn(f)
        warning('GAIL:cubSobol_SI_fo_g:fnotfcn',...
            'The given input f was not a function. Example for f(x)=x^2:')
        f = @(x) x.^2;
        out_param.f=f;
        hyperbox = default.hyperbox;
    else
        out_param.f=f;
        hyperbox = varargin{2};
        if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(size(hyperbox,2)<370)
            warning('GAIL:cubSobol_SI_fo_g:hyperbox_error1',...
                'The hyperbox must be a real matrix of size 2xd where d can not be greater than 370. Example for f(x)=x^2:')
            f = @(x) x.^2;
            out_param.f=f;
            hyperbox = default.hyperbox;
        end
    end
end

validvarargin=numel(varargin)>2;
if validvarargin
    in3=varargin(3:end);
    for j=1:numel(varargin)-2
    validvarargin=validvarargin && (isnumeric(in3{j}) ...
        || ischar(in3{j}) || isstruct(in3{j}) || gail.isfcn(in3{j}));
    end
    if ~validvarargin
        warning('GAIL:cubSobol_SI_fo_g:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in3=varargin{3};
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
    if isnumeric(in3) || ischar(in3)
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
        if isstruct(in3) %parse input structure
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
    parse(p,f,hyperbox,varargin{3:end})
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
if ~isnumeric(hyperbox) || ~(size(hyperbox,1)==2) || ~(out_param.d<370)
    warning('GAIL:cubSobol_SI_fo_g:hyperbox_error2',...
        'The hyperbox must be a real matrix of size 2 x d where d can not be greater than 370. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    hyperbox = default.hyperbox;
end

% Force measure to be uniform or normal only
if ~(strcmp(out_param.measure,'uniform') || strcmp(out_param.measure,'normal') )
    warning('GAIL:cubSobol_SI_fo_g:notmeasure',['The measure can only be uniform or normal.' ...
            ' Using default measure ' num2str(default.measure)])
    out_param.measure = default.measure;
end

% Force absolute tolerance greater than 0
if (out_param.abstol < 0 )
    warning('GAIL:cubSobol_SI_fo_g:abstolnonpos',['Absolute tolerance cannot be negative.' ...
            ' Using default absolute tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force relative tolerance greater than 0 and smaller than 1
if (out_param.reltol < 0) || (out_param.reltol > 1)
    warning('GAIL:cubSobol_SI_fo_g:reltolnonunit',['Relative tolerance should be chosen in [0,1].' ...
            ' Using default relative tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('GAIL:cubSobol_SI_fo_g:lowmmin',['The minimum starting exponent ' ...
            'should be an integer greater than 0 and smaller or equal than the maxium.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force mmin to be integer greater than r_lag (so that l_star=mmin-r_lag>=0)
if out_param.mmin < r_lag
    warning('GAIL:cubSobol_SI_fo_g:lowmminrlag',['The minimum starting exponent ' ...
            'should be at least ' num2str(r_lag) '.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 54
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=53)
    warning('GAIL:cubSobol_SI_fo_g:wrongmmax',['The maximum exponent for the budget should be an integer biger than mmin and smaller than 54.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('GAIL:cubSobol_SI_fo_g:fudgenonpos',['The fudge factor should be a positive function.' ...
            ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end

% Force toltype to be max or comb
if ~(strcmp(out_param.toltype,'max') || strcmp(out_param.toltype,'comb') )
    warning('GAIL:cubSobol_SI_fo_g:nottoltype',['The error type can only be max or comb.' ...
            ' Using default toltype ' num2str(default.toltype)])
    out_param.toltype = default.toltype;
end

% Force theta to be in [0,1]
if (out_param.theta < 0) || (out_param.theta > 1)
    warning('GAIL:cubSobol_SI_fo_g:thetanonunit',['Theta should be chosen in [0,1].' ...
            ' Using default theta ' num2str(default.theta)])
    out_param.theta = default.theta;
end

% Checking on pure absolute/relative error
if (out_param.abstol==0) && (out_param.reltol==0)
    warning('GAIL:cubSobol_SI_fo_g:tolzeros',['Absolute and relative error tolerances can not be simultaniusly 0.' ...
            ' Using default absolute tolerance ' num2str(default.abstol) ' and relative tolerance ' num2str(default.reltol)])
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
end
if (strcmp(out_param.toltype,'comb')) && (out_param.theta==1) && (out_param.abstol==0)
    warning('GAIL:cubSobol_SI_fo_g:abstolzero',['When choosing toltype comb, if theta=1 then abstol>0.' ...
            ' Using default absolute tolerance ' num2str(default.abstol) ])
    out_param.abstol = default.abstol;
end
if (strcmp(out_param.toltype,'comb')) && (out_param.theta==0) && (out_param.reltol==0)
    warning('GAIL:cubSobol_SI_fo_g:reltolzero',['When choosing toltype comb, if theta=0 then reltol>0.' ...
            ' Using default relative tolerance ' num2str(default.reltol) ])
    out_param.reltol = default.reltol;
end

% Checking on the hyperbox given the measure
if (strcmp(out_param.measure,'uniform')) && ~all(all(isfinite(hyperbox)))
    warning('GAIL:cubSobol_SI_fo_g:hyperboxnotfinite',['If uniform measure, hyperbox must be of finite volume.' ...
            ' Using default hyperbox:'])
    disp([zeros(1,out_param.d);ones(1,out_param.d)])
    hyperbox = [zeros(1,out_param.d);ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (sum(sum(isfinite(hyperbox)))>0)
    warning('GAIL:cubSobol_SI_fo_g:hyperboxfinite',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end
if (strcmp(out_param.measure,'normal')) && (any(hyperbox(1,:)==hyperbox(2,:)) || any(hyperbox(1,:)>hyperbox(2,:)))
    warning('GAIL:cubSobol_SI_fo_g:hyperboxnormalwrong',['If normal measure, hyperbox must be defined as (-Inf,Inf)^d.' ...
            ' Using default hyperbox:'])
    disp([-inf*ones(1,out_param.d);inf*ones(1,out_param.d)])
    hyperbox = [-inf*ones(1,out_param.d);inf*ones(1,out_param.d)];
end

end


function y = estimate_min(down, up)
    if down(4)*up(4) <= 0
        y = min(max(max(down(1)-max(down(2).^2, up(2).^2)...
        ,0)./max(up(3), eps), 0),1);
    else
        y = min(max(max(down(1)-max(down(2).^2, up(2).^2)...
        ,0)./max(up(3)-min(down(4).^2, up(4).^2), eps), 0),1);
    end
end

function y = estimate_max(down, up)
    if down(2)*up(2) <= 0
        y = min(max(max(up(1)...
    ,0)./max(down(3)-max(down(4).^2, up(4).^2), eps), 0),1);
    else
        y = min(max(max(up(1)-min(down(2).^2, up(2).^2)...
    ,0)./max(down(3)-max(down(4).^2, up(4).^2), eps), 0),1);
    end
end