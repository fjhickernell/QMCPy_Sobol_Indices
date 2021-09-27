clearvars

abstol = 2*1e-2;
reltol = 0e-2; % Pure absolute tolerance
mmin = 10;
mmax = 20; % I adjust that not to run out of memory. It can go up to 54. Type help cubSobol_SI_g for more information.
wronga = 0; % Wrong estimates with all indices method
wrong = 0; % Wrong estimates using the one by one algo
s = .0; %time in seconds used delay each function evaluation

fudge = @(m,d) 5*2.^-(m);

%% Ishigami
disp('Running Ishigami example ...')
d = 3;
f = @(x) sin(x(:,1)).*(1+1/10*(x(:,3).^4))+7*sin(x(:,2)).^2 +  pause_t(s);
hyperbox = pi*[-ones(1,d) ; ones(1,d)];

R = [.3139051827, .4424111333, 0.; .5575888667, .4424111326, .2436836832];
[q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
SI = q;
SI_n = out_param.n;
SI_t = out_param.time;
disp(SI)
% disp(out_param.bound_err)
disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > max(abstol,R*reltol))
% disp(out_param.exitflag)
% wronga = wronga + sum(sum(abs(R-SI) > max(abstol,R*reltol)));
% disp(['Time:' num2str(SI_t)])

% disp('Running Ishigami example ...')
% exitflag = [];
% exitflagt = [];
% for j = 1:d
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflag = [exitflag out_param.exitflag];
%     SI(1,j) = q;
%     SI_n(1,j) = out_param.n;
%     SI_t(1,j) = out_param.time;
%     SI_s(1,j) = out_param.small;
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflagt = [exitflagt out_param.exitflag];
%     SI(2,j) = q;
%     SI_n(2,j) = out_param.n;
%     SI_t(2,j) = out_param.time;
%     SI_s(2,j) = out_param.small;
% end
% disp(SI)
% disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > abstol)
% wrong = wrong + sum(sum(abs(R-SI) > abstol));
% if any(sum(exitflag) > 0) || any(sum(exitflagt) > 0)
%     warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
% end
% disp(['Time:' num2str(sum(sum(SI_t)))])


%% Bratley et al.
disp('Running Bratley et al. example ...')
d = 6;
f = @(x) sum(bsxfun(@times, cumprod(x,2), (-1).^(1:d)),2) +  pause_t(s);
hyperbox = [-zeros(1,d) ; ones(1,d)];

R = [.6528636616, .1791303924, 0.3701041165e-1, 0.1332374820e-1, 0.1480416466e-2, 0.1480416466e-2; .7396477462, .2659144770, 0.764211693e-1, 0.343115454e-1, 0.62384628e-2, 0.62384628e-2];
[q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
SI = q;
SI_n = out_param.n;
SI_t = out_param.time;
disp(SI)
% disp(out_param.bound_err)
disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > max(abstol,R*reltol))
% wronga = wronga + sum(sum(abs(R-SI) > max(abstol,R*reltol)));
% disp(['Time:' num2str(SI_t)])

% disp('Running Bratley et al. example ...')
% exitflag = [];
% exitflagt = [];
% for j = 1:d
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflag = [exitflag out_param.exitflag];
%     SI(1,j) = q;
%     SI_n(1,j) = out_param.n;
%     SI_t(1,j) = out_param.time;
%     SI_s(1,j) = out_param.small;
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflagt = [exitflagt out_param.exitflag];
%     SI(2,j) = q;
%     SI_n(2,j) = out_param.n;
%     SI_t(2,j) = out_param.time;
%     SI_s(2,j) = out_param.small;
% end
% disp(SI)
% disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > abstol)
% wrong = wrong + sum(sum(abs(R-SI) > abstol));
% if any(sum(exitflag) > 0) || any(sum(exitflagt) > 0)
%     warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
% end
% disp(['Time:' num2str(sum(sum(SI_t)))])



%% Sobol' g-function
disp('Running Sobol g-function example ...')
d = 6;
a = [0 1/2 3 9 99 99];
f = @(x) prod(bsxfun(@(x,a) (abs(4*x-2)+a)./(1+a), x , a),2) +  pause_t(s);
hyperbox = [zeros(1,d) ; ones(1,d)];

R = [.5867811893, .2607916397, 0.3667382378e-1, 0.5867811312e-2, 0.5867753221e-4, 0.5867753221e-4; .6900858920, .3561733634, 0.563335432e-1, 0.91705767e-2, 0.920079e-4, 0.920079e-4];
[q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
SI = q;
SI_n = out_param.n;
SI_t = out_param.time;
disp(SI)
% disp(out_param.bound_err)
disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > max(abstol,R*reltol))
% wronga = wronga + sum(sum(abs(R-SI) > max(abstol,R*reltol)));
% disp(['Time:' num2str(SI_t)])

% disp('Running Sobol g-function example ...')
% exitflag = [];
% exitflagt = [];
% R = [.5867811893, .2607916397, 0.3667382378e-1, 0.5867811312e-2, 0.5867753221e-4, 0.5867753221e-4; .6900858920, .3561733634, 0.563335432e-1, 0.91705767e-2, 0.920079e-4, 0.920079e-4];
% for j = 1:d
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflag = [exitflag out_param.exitflag];
%     SI(1,j) = q;
%     SI_n(1,j) = out_param.n;
%     SI_t(1,j) = out_param.time;
%     SI_s(1,j) = out_param.small;
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflagt = [exitflagt out_param.exitflag];
%     SI(2,j) = q;
%     SI_n(2,j) = out_param.n;
%     SI_t(2,j) = out_param.time;
%     SI_s(2,j) = out_param.small;
% end
% disp(SI)
% disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > abstol)
% wrong = wrong + sum(sum(abs(R-SI) > abstol));
% if any(sum(exitflag) > 0) || any(sum(exitflagt) > 0)
%     warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
% end
% disp(['Time:' num2str(sum(sum(SI_t)))])


%% Morokoff and Caflish
disp('Running Morokoff and Caflish example ...')
d = 6;
f = @(x) (1+1/d)^d*prod(x,2).^(1/d) +  pause_t(s);
hyperbox = [zeros(1,d) ; ones(1,d)];

R = [.1581948744, .1581948744, .1581948744, .1581948744, .1581948744, .1581948744; .1753745708, .1753745708, .1753745708, .1753745708, .1753745708, .1753745708];
[q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
SI = q;
SI_n = out_param.n;
SI_t = out_param.time;
disp(SI)
% disp(out_param.bound_err)
disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > max(abstol,R*reltol))
% wronga = wronga + sum(sum(abs(R-SI) > max(abstol,R*reltol)));
% disp(['Time:' num2str(SI_t)])
% disp(['Wrong estimates with all:' num2str(wronga)]);

% disp('Running Morokoff and Caflish example ...')
% exitflag = [];
% exitflagt = [];
% for j = 1:d
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflag = [exitflag out_param.exitflag];
%     SI(1,j) = q;
%     SI_n(1,j) = out_param.n;
%     SI_t(1,j) = out_param.time;
%     SI_s(1,j) = out_param.small;
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflagt = [exitflagt out_param.exitflag];
%     SI(2,j) = q;
%     SI_n(2,j) = out_param.n;
%     SI_t(2,j) = out_param.time;
%     SI_s(2,j) = out_param.small;
% end
% disp(SI)
% disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > abstol)
% wrong = wrong + sum(sum(abs(R-SI) > abstol));
% if any(sum(exitflag) > 0) || any(sum(exitflagt) > 0)
%     warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
% end
% disp(['Time:' num2str(sum(sum(SI_t)))])
% disp(['Wrong estimates one by one:' num2str(wrong)]);

%% Example Laurent
disp('Running Laurents example ...')
d = 6;
f = @(x) (x(:,3).^3).*exp(x(:,4)).*(x(:,5).^2).*(x(:,2).^3.*x(:,3)+sin(x(:,4)+1).*x(:,6))./(1+x(:,1)).^(1/2);
hyperbox = [zeros(1,d) ; ones(1,d)];

R = [0.1949202354e-2, 0.2167739767e-1, .2765283730, 0.1658947475e-1, .1561007328, 0.3242382766e-1 0.118205783e-1, .1184436834, .7035364275, 0.934641572e-1, .5311670735, .1464118928];
[q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
SI = q;
SI_n = out_param.n;
SI_t = out_param.time;
disp(SI)
% disp(out_param.bound_err)
disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > max(abstol,R*reltol))
% wronga = wronga + sum(sum(abs(R-SI) > max(abstol,R*reltol)));
% disp(['Time:' num2str(SI_t)])
% disp(['Wrong estimates with all:' num2str(wronga)]);

% disp('Running Morokoff and Caflish example ...')
% exitflag = [];
% exitflagt = [];
% for j = 1:d
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflag = [exitflag out_param.exitflag];
%     SI(1,j) = q;
%     SI_n(1,j) = out_param.n;
%     SI_t(1,j) = out_param.time;
%     SI_s(1,j) = out_param.small;
%     [q,app_int,out_param] = cubSobol_SI_g(f,hyperbox,-j,'measure','uniform','abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',fudge);
%     exitflagt = [exitflagt out_param.exitflag];
%     SI(2,j) = q;
%     SI_n(2,j) = out_param.n;
%     SI_t(2,j) = out_param.time;
%     SI_s(2,j) = out_param.small;
% end
% disp(SI)
% disp(SI_n)
% disp(SI_s)
% disp(abs(R-SI) > abstol)
% wrong = wrong + sum(sum(abs(R-SI) > abstol));
% if any(sum(exitflag) > 0) || any(sum(exitflagt) > 0)
%     warning('Results cannot be guaranteed due to reaching maximum budget or failing necessary conditions.')
% end
% disp(['Time:' num2str(sum(sum(SI_t)))])
% disp(['Wrong estimates one by one:' num2str(wrong)]);