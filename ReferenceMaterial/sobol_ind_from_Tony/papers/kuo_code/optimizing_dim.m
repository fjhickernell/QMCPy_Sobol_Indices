s = 25; % Number of matrices for which we generate the Sobol' seq generator
m = 20; % Number of digits used to generate the seq up to 2^m points

weight_matrix = (s-1)*ones(2*s,2*s);
weight_matrix(1:s,1:s) = s;
weight_matrix(s+1:end,s+1:end) = s-2;
weight_matrix(1:2*s+1:4*s*s) = 0;
% for i = 1:2*s
%     for j = 1:2*s
%         if j >= i
%             weight_matrix(i,j) = 0;
%         end
%     end
% end

table_t_values_file = sprintf('tval_table_s%d_m%d.mat',2*s,m);
if exist(table_t_values_file)
    load(table_t_values_file)
else
    error('No t-values table found to optimize.')
end

initial_value = sum(sum(t_values_table.*weight_matrix));
fname1 = sprintf('better_permutations_s%d_m%d.mat',s,m);
fname2 = sprintf('weight_values_s%d_m%d.mat',s,m);
if exist(fname1)
    load(fname1)
end
if exist(fname2)
    load(fname2)
end
better_values = [];
new_permutations = [];
fraction = 1/2; % fraction of s dimensions that we keep their order at the beginning
tic
parfor i = 1:20000000
    p = [1:floor(s*fraction) floor(s*fraction)+randperm(s+ceil(s*(1-fraction)))];
    new_value = sum(sum(t_values_table(p,p).*weight_matrix));
    if initial_value > new_value;
        disp(new_value)
        better_values = [better_values; new_value];
        new_permutations = [new_permutations; p];
    end
end
toc

values = [values; better_values];
permutations = [permutations; new_permutations];
[~,I] = sort(values); % we only keep the best n weighted values
n = min(size(values,1),5); % we keep the best 5 values or less if we have less
values = values(I(1:n));
permutations = permutations(I(1:n),:);
save(fname1, 'permutations')
save(fname2, 'values')