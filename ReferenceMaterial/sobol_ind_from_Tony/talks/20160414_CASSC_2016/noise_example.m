d = 3;
f = @(u) exp(1.5*u(:,1)+u(:,2)) + cos(2*pi*u(:,3))
hyperbox = [-zeros(1,d); ones(1,d)];

mmin = 15;
mmax = 20;
abstol = 1e-4;
reltol = 0e-3;
measure = 'uniform';

[q_mean, out] = cubSobol_g(f, hyperbox, measure, abstol, reltol, 'mmin', mmin, 'mmax', mmax);
disp(q_mean)
[q, out] = cubSobol_SI_all_g(f, hyperbox, measure, abstol, reltol, 'mmin', mmin, 'mmax', mmax);
disp(q)
[q_var, out] = cubSobol_g(@(u) f(u).^2 - q_mean^2, hyperbox, measure, abstol, reltol, 'mmin', mmin, 'mmax', mmax);
disp(q_var)