clearvars

abstol = 5*1e-3;
pre = 2; % Precision to display according to 1*1e-3 (3 digits)
reltol = 0e-2; % Pure absolute tolerance
mmin = 9;
mmax = 24; % I adjust that not to run out of memory. It can go up to 54. Type help cubSobol_SI_g for more information.
threshold_small = 0.1; % Below which sizes do we use correl2 estimator
s = .0; %time in seconds used delay each function evaluation
samples = 100;

fudge = @(m,d) 10*2.^-(1.*m);

% %% Bratley et al.
% disp('Running Bratley et al. example ...')
% d = 6;
% f = @(x) sum(bsxfun(@times, cumprod(x,2), (-1).^(1:d)),2) +  pause_t(s);
% hyperbox = [-zeros(1,d) ; ones(1,d)];
% 
% R = [.6528636616, .1791303924, 0.3701041165e-1, 0.1332374820e-1, 0.1480416466e-2, 0.1480416466e-2; .7396477462, .2659144770, 0.764211693e-1, 0.343115454e-1, 0.62384628e-2, 0.62384628e-2];
% SI = zeros(2, d);
% SI_small = zeros(2, d);
% SI_n = SI;
% SI_n_small = SI_small;
% error = SI;
% error_small = SI_small;
% reliability = SI;
% reliability_small = SI_small;
% error_all = [];
% error_all_small = [];
% SI_n_print = [];
% SI_n_print_small = [];
% small = [];
% for k = 1:samples
%     sobstr = sobolset(3*d); %generate a Sobol' sequence 3*d to consider the changing to the estimator for some smaller size indices
%     sobstr = scramble(sobstr,'MatousekAffineOwen'); %scramble it
%     [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,sobstr,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', 0);
%     SI = SI + q;
%     SI_n_print = [SI_n_print; out_param.n];
%     SI_n = SI_n + out_param.n;
%     error = error + abs(R-q);
%     error_all = [error_all; abs(R-q)];
%     reliability = reliability + 1*(abs(R-q) > max(abstol,R*reltol));
%     [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,sobstr,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', threshold_small);
%     SI_small = SI_small + q;
%     SI_n_print_small = [SI_n_print_small; out_param.n];
%     SI_n_small = SI_n_small + out_param.n;
%     error_small = error_small + abs(R-q);
%     error_all_small = [error_all_small; abs(R-q)];
%     reliability_small = reliability_small + 1*(abs(R-q) > max(abstol,R*reltol));
%     small = [small; out_param.small];
% end
% SI = SI/samples; SI_n = SI_n/samples; reliability = 1 - reliability/samples; error = error/samples;
% SI_small = SI_small/samples; SI_n_small = SI_n_small/samples; reliability_small = 1 - reliability_small/samples; error_small = error_small/samples;
% 
% csvwrite('bratley_variantAa_reliability.csv', reliability)
% csvwrite('bratley_variantAa_errors.csv', error_all)
% csvwrite('bratley_variantAa_nvalues.csv', SI_n_print)
% csvwrite('bratley_variantAb_reliability.csv', reliability_small)
% csvwrite('bratley_variantAb_errors.csv', error_all_small)
% csvwrite('bratley_variantAb_nvalues.csv', SI_n_print_small)
% csvwrite('bratley_variantAb_small_estimates.csv', small)
% 
% 
% %% Sobol' g-function
% disp('Running Sobol g-function example ...')
% d = 6;
% a = [0 1/2 3 9 99 99];
% f = @(x) prod(bsxfun(@(x,a) (abs(3*x-2)+a)./(1+a), x , a),2) +  pause_t(s);
% hyperbox = [zeros(1,d) ; ones(1,d)];
% 
% % R = [.5867811893, .2607916397, 0.3667382378e-1, 0.5867811312e-2, 0.5867753221e-4, 0.5867753221e-4; .6900858920, .3561733634, 0.563335432e-1, 0.91705767e-2, 0.920079e-4, 0.920079e-4]; % If abs(4x-2)
% R = [.6042800971, .2360469148, 0.2855765831e-1, 0.4339846615e-2, 0.4210460219e-4, 0.4210460219e-4; .7251945137, .3480933667, 0.483463013e-1, 0.74762327e-2, 0.727601e-4, 0.727601e-4]; % If abs(3x-2)
% SI = zeros(2, d);
% SI_small = zeros(2, d);
% SI_n = SI;
% SI_n_small = SI_small;
% error = SI;
% error_small = SI_small;
% reliability = SI;
% reliability_small = SI_small;
% error_all = [];
% error_all_small = [];
% SI_n_print = [];
% SI_n_print_small = [];
% small = [];
% for k = 1:samples
%     disp(k)
%     sobstr = sobolset(3*d); %generate a Sobol' sequence 3*d to consider the changing to the estimator for some smaller size indices
%     sobstr = scramble(sobstr,'MatousekAffineOwen'); %scramble it
%     [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,sobstr,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', 0);
%     SI = SI + q;
%     SI_n_print = [SI_n_print; out_param.n];
%     SI_n = SI_n + out_param.n;
%     error = error + abs(R-q);
%     error_all = [error_all; abs(R-q)];
%     reliability = reliability + 1*(abs(R-q) > max(abstol,R*reltol));
%     [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,sobstr,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', threshold_small);
%     SI_small = SI_small + q;
%     SI_n_print_small = [SI_n_print_small; out_param.n];
%     SI_n_small = SI_n_small + out_param.n;
%     error_small = error_small + abs(R-q);
%     error_all_small = [error_all_small; abs(R-q)];
%     reliability_small = reliability_small + 1*(abs(R-q) > max(abstol,R*reltol));
%     small = [small; out_param.small];
% end
% SI = SI/samples; SI_n = SI_n/samples; reliability = 1 - reliability/samples; error = error/samples;
% SI_small = SI_small/samples; SI_n_small = SI_n_small/samples; reliability_small = 1 - reliability_small/samples; error_small = error_small/samples;
% 
% csvwrite('sobolg_variantAa_reliability.csv', reliability)
% csvwrite('sobolg_variantAa_errors.csv', error_all)
% csvwrite('sobolg_variantAa_nvalues.csv', SI_n_print)
% csvwrite('sobolg_variantAb_reliability.csv', reliability_small)
% csvwrite('sobolg_variantAb_errors.csv', error_all_small)
% csvwrite('sobolg_variantAb_nvalues.csv', SI_n_print_small)
% csvwrite('sobolg_variantAb_small_estimates.csv', small)


%% Wing weight function
disp('Running wing weight function example ...')
d = 15;
f = @(x) 0.036*x(:, 1).^0.758.*x(:, 2).^0.0035.*(x(:, 3)./ ...
    cos(x(:, 4)).^2).^0.6.*x(:, 5).^0.006.*x(:, 6).^0.04 ...
    .*(100*x(:, 7)./cos(x(:, 4))).^-0.3.*(x(:, 8).*x(:, 9)).^0.49 + ...
    x(:, 1).*x(:, 10);
hyperbox = [150 220 6  -1/18*pi 16 0.5 0.08 2.5 1700 0.025 zeros(1, d-10);
            200 300 10  1/18*pi 45 1   0.18 6   2500 0.08 ones(1, d-10)];

SI = zeros(2, d);
SI_estimates = [];
SI_estimates_small = [];
SI_small = zeros(2, d);
SI_n = SI;
SI_n_small = SI_small;
SI_n_print = [];
SI_n_print_small = [];
small = [];
for k = 1:samples
    sobstr = sobolset(3*d); %generate a Sobol' sequence 3*d to consider the changing to the estimator for some smaller size indices
    sobstr = scramble(sobstr,'MatousekAffineOwen'); %scramble it
    [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,sobstr,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', 0);
    SI = SI + q;
    SI_estimates = [SI_estimates; q];
    SI_n_print = [SI_n_print; out_param.n];
    SI_n = SI_n + out_param.n;
    [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,sobstr,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', threshold_small);
    SI_small = SI_small + q;
    SI_estimates_small = [SI_estimates_small; q];
    SI_n_print_small = [SI_n_print_small; out_param.n];
    SI_n_small = SI_n_small + out_param.n;
    small = [small; out_param.small];
end
SI = SI/samples; SI_n = SI_n/samples;
SI_small = SI_small/samples; SI_n_small = SI_n_small/samples;

csvwrite('wingweight_variantAa_nvalues.csv', SI_n_print)
csvwrite('wingweight_variantAa_SIvalues.csv', SI_estimates)
csvwrite('wingweight_variantAb_SIvalues.csv', SI_estimates_small)
csvwrite('wingweight_variantAb_nvalues.csv', SI_n_print_small)
csvwrite('wingweight_variantAb_small_estimates.csv', small)
