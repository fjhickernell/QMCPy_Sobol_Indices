% %% Construction of C generating matrices
% s = 28; % Number of matrices for which we generate the Sobol' seq generator
% m = 12; % Number of digits used to generate the seq up to 2^m points
% entry_file = 'new-joe-kuo-6.21201'; % joe-kuo-old.1111 new-joe-kuo-6.21201
% A = txt2mat(entry_file);
% 
% %% We take the generator values to the struct C and generate matrices Cs
% for k = 1:size(A,1)
%     C(k + 1).d = A(k,1); % Dimension
%     C(k + 1).s = A(k,2); % Polynomial degree
%     C(k + 1).a = A(k,3); % Polynomial coefficients assuming a = a_d-1 * 2^0 + a_d-2 * 2^1 + ... + a_1 2^d-2
%     aa = find(isnan(A(k,:))==1);
%     if ~isempty(aa)
%         C(k + 1).m = A(k,4:aa(1) - 1); % Initial directional numbers m_i
%     else
%         C(k + 1).m = A(k,4:end);
%     end
% end
% 
% C(1).C = eye(m); C(1).d = 1; C(1).s = 0; C(1).a = []; C(1).m = 1;
% C(2).C = C(1).C;
% % I generate C(2).C manually
% for i = 2:m
%     C(2).m = [C(2).m bitxor(2*C(2).m(end),C(2).m(end))];
%     aa = (dec2bin(C(2).m(i))- '0')';
%     empty = i - size(aa,1);
%     C(2).C(1:i,i) = [zeros(empty,1); aa];
% end
% 
% % Generate C(j).C matrices for j = 3:s
% for j = 3:s % For each dimension
%     C(j).C = C(1).C;
%     for i = 2:m % for each column of C(j).C
%         if i > size(C(j).m,2)
%             vec_a = dec2bin(C(j).a)- '0';
%             if size(vec_a,2) < C(j).s - 1
%                 vec_a = [zeros(1,C(j).s-size(vec_a,2)-1), vec_a];
%             end
%             vec_a = [(vec_a), 1, 1].*[2.^(1:C(j).s), 1]; % Matlab mistake: [fliplr(vec_a), 1, 1].*[2.^(1:C(j).s), 1]
%             vec_a = vec_a.*[C(j).m(size(C(j).m,2):-1:size(C(j).m,2)-C(j).s+1), C(j).m(size(C(j).m,2)-C(j).s+1)];
%             C(j).m = [C(j).m, bitxor_toni(vec_a)];
%         end
%         aa = (dec2bin(C(j).m(i))- '0')';
%         empty = i - size(aa,1);
%         C(j).C(1:i,i) = [zeros(empty,1);aa];
%     end
% end
% 
% 
%% We study the t-values of the sequence
d = 4;
dim_set = [4     6     1     2     5     3     8     7];%randperm(2*d); %1:2*d   [4     5     8     6     3     7     1     2]
%[4    17    12    18     1     2    15     9    11    14    16     3     8      7    13    10     5    19     6     20]
% [10     8     3     9     6     4     5     2     1     7]
% dim_set = [randperm(d), d + randperm(d)];

clear average_t_values max_t_values
for k = 1:d
    clear SEQ
    U = mod(C(dim_set(d+k)).C\C(dim_set(k)).C, 2);
    parfor i = 1:2*d - 1
        if i < d + 1
            SEQ(i).C = C(dim_set(i)).C;
        elseif i < d + k
            SEQ(i).C = mod(C(dim_set(i)).C*U, 2);
        else
            SEQ(i).C = mod(C(dim_set(i + 1)).C*U, 2);
        end
    end

    %% Evaluating t-values and creating the t-values table for SEQ
    t_values_table = zeros(2*d - 1, 2*d - 1);
    for dim = 1:2*d - 1
        for dim2 = 1:2*d - 1
            if dim <= dim2;
                t_values_table(dim,dim2) = 0;
            else
                t_values_table(dim,dim2) = t_value(SEQ(dim).C, SEQ(dim2).C, m);
            end
        end
%         disp(t_values_table(dim,1:dim-1))
    end
    average_t_values(k) = sum(sum(t_values_table))/((2*d - 1)*(2*d - 2)/2);
    max_t_values(k) = max(max(t_values_table));
end
disp([dim_set])
disp([num2str(average_t_values) ' -> ' num2str(mean(average_t_values))])
disp([num2str(max_t_values) ' -> ' num2str(mean(max_t_values))])