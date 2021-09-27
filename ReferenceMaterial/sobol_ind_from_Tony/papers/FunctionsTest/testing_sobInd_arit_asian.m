abstol = 5*1e-3;
pre = 3; % Precision to display according to 1*1e-3 (3 digits)
reltol = 0e-2; % Pure absolute tolerance
mmin = 9;
mmax = 24; % I adjust that not to run out of memory. It can go up to 54. Type help cubSobol_SI_g for more information.
samples = 100;
fudge = @(m,d) 10*2.^-(1.*m);
d = 6;

inp.timeDim.timeVector = 1/52:1/52:d*1/52; %weekly monitoring for 1 year
inp.assetParam.initPrice = 36; %initial stock price
inp.assetParam.interest = 0.06; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 40; %strike price
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'PCA';
EuroCall = optPrice(inp); %construct an optPrice object 

opt = optPayoff(EuroCall); %make a copy
opt.payoffParam = struct('optType',{{'amean'}}, 'putCallType',{{'call'}}); 

f =@(x) genOptPayoffs(opt,x);
hyperbox = [zeros(1,d) ; ones(1,d)];

%% Testing without small indices estimator
SI = zeros(2, d);
SI_n = SI;
tic
for k = 1:samples
    [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', 0);
    SI = SI + q;
    SI_n = SI_n + out_param.n;
end
toc
SI = SI/samples; SI_n = SI_n/samples;
round(SI, pre, 'significant')
round(SI_n)

%% Testing with small indices estimator
SI = zeros(2, d);
SI_n = SI;
tic
for k = 1:samples
    [q,app_int,out_param] = cubSobol_SI_all_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d), 'threshold_small', 0.1);
    SI = SI + q;
    SI_n = SI_n + out_param.n;
end
toc
SI = SI/samples; SI_n = SI_n/samples;
round(SI, pre, 'significant')
round(SI_n)

%% Testing replicated
SI = zeros(1, d);
SI_n = SI;
tic
for k = 1:samples
    [q,app_int,out_param] = cubSobol_SI_fo_g(f,hyperbox,'abstol',abstol,'reltol',reltol,'mmin',mmin,'mmax',mmax,'fudge',@(m) fudge(m,d));
    SI = SI + q;
    SI_n = SI_n + out_param.n;
end
toc
SI = SI/samples; SI_n = SI_n/samples;
round(SI, pre, 'significant')
round(SI_n)