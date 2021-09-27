clearvars

abstol = 5*1e-3;
reltol = 0e-2; % Pure absolute tolerance
mmin = 9;
mmax = 24; % I adjust that not to run out of memory. It can go up to 54. Type help cubSobol_SI_g for more information.
wronga = 0; % Wrong estimates with all indices method
wrong = 0; % Wrong estimates using the one by one algo
s = .0; %time in seconds used delay each function evaluation
samples = 100;

fudge = @(m,d) 10*2.^-(1.*m);

% %% Bratley et al.
% disp('Running Bratley et al. example ...')
% d = 6;
% f = @(x) sum(bsxfun(@times, cumprod(x,2), (-1).^(1:d)),2) +  pause_t(s);
% hyperbox = [-zeros(1,d) ; ones(1,d)];
% 
% R = [.6528636616, .1791303924, 0.3701041165e-1, 0.1332374820e-1, 0.1480416466e-2, 0.1480416466e-2];
% SI = zeros(1, d);
% SI_n = SI;
% error = SI;
% reliability = SI;
% error_all = [];
% SI_n_print = [];
% for k = 1:samples
%     [q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
%     SI = SI + q;
%     SI_n_print = [SI_n_print; out_param.n];
%     SI_n = SI_n + out_param.n;
%     error = error + abs(R-q);
%     error_all = [error_all; abs(R-q)];
%     reliability = reliability + 1*(abs(R-q) > max(abstol,R*reltol));
% end
% SI = SI/samples; SI_n = SI_n/samples; reliability = 1 - reliability/samples; error = error/samples;
% % round(SI, pre, 'significant')
% % round(SI_n)
% % round(error, pre, 'significant')
% % round(reliability*100)
% 
% csvwrite('bratley_variantB_replicated_reliability.csv', reliability)
% csvwrite('bratley_variantB_replicated_errors.csv', error_all)
% csvwrite('bratley_variantB_replicated_nvalues.csv', SI_n_print)
% 
% 

% %% Sobol' g-function
% disp('Running Sobol g-function example ...')
% d = 6;
% a = [0 1/2 3 9 99 99];
% f = @(x) prod(bsxfun(@(x,a) (abs(3*x-2)+a)./(1+a), x , a),2) +  pause_t(s);
% hyperbox = [zeros(1,d) ; ones(1,d)];
% 
% % R = [.5867811893, .2607916397, 0.3667382378e-1, 0.5867811312e-2, 0.5867753221e-4, 0.5867753221e-4]; % If abs(4x-2)
% R = [.6042800971, .2360469148, 0.2855765831e-1, 0.4339846615e-2, 0.4210460219e-4, 0.4210460219e-4]; % If abs(3x-2)
% SI = zeros(1, d);
% SI_n = SI;
% error = SI;
% reliability = SI;
% error_all = [];
% SI_n_print = [];
% for k = 1:samples
%     [q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
%     SI = SI + q;
%     SI_n_print = [SI_n_print; out_param.n];
%     SI_n = SI_n + out_param.n;
%     error = error + abs(R-q);
%     error_all = [error_all; abs(R-q)];
%     reliability = reliability + 1*(abs(R-q) > max(abstol,R*reltol));
% end
% SI = SI/samples; SI_n = SI_n/samples; reliability = 1 - reliability/samples; error = error/samples;
% % round(SI, pre,'significant')
% % round(SI_n)
% % round(error, pre, 'significant')
% % round(reliability*100)
% 
% csvwrite('sobolg_variantB_replicated_reliability.csv', reliability)
% csvwrite('sobolg_variantB_replicated_errors.csv', error_all)
% csvwrite('sobolg_variantB_replicated_nvalues.csv', SI_n_print)
% 
% 

%% Wing weight function
disp('Running wing weight function example ...')
d = 15;
f = @(x) 0.036*x(:, 1).^0.758.*x(:, 2).^0.0035.*(x(:, 3)./ ...
    cos(x(:, 4)).^2).^0.6.*x(:, 5).^0.006.*x(:, 6).^0.04 ...
    .*(100*x(:, 7)./cos(x(:, 4))).^-0.3.*(x(:, 8).*x(:, 9)).^0.49 + ...
    x(:, 1).*x(:, 10);
hyperbox = [150 220 6  -1/18*pi 16 0.5 0.08 2.5 1700 0.025 zeros(1, d-10);
            200 300 10  1/18*pi 45 1   0.18 6   2500 0.08 ones(1, d-10)];

SI = zeros(1, d);
SI_n = SI;
SI_n_print = [];
SI_estimates = [];
for k = 1:samples
    [q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
    SI = SI + q;
    SI_estimates = [SI_estimates; q];
    SI_n_print = [SI_n_print; out_param.n];
    SI_n = SI_n + out_param.n;
end
SI = SI/samples; SI_n = SI_n/samples;

csvwrite('wingweight_variantB_replicated_nvalues.csv', SI_n_print)
csvwrite('wingweight_variantB_replicated_SIvalues.csv', SI_estimates)