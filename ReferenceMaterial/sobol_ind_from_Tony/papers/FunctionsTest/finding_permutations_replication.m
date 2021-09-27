tic
n = 2^17;
d = 100;
sob = sobolset(d);

x = net(sob,n);
perm = zeros(n,d);
perm_inv = perm;
for i = 1:d
    [~, perm(:,i)] = sort(x(:,i));
    [~, perm_inv(:,i)] = sort(perm(:,i));
end
toc
tic
save('perm.mat','perm','perm_inv','-v7.3')
toc
