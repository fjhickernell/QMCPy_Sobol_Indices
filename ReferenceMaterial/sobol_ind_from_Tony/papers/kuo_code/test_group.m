clearvars
m = 7;
% m:2   max group order:2   group size:2       2^1
% m:3   max group order:4   group size:8       2^3
% m:4   max group order:4   group size:64      2^6
% m:5   max group order:8   group size:1024    2^10
% m:6   max group order:8   group size:32768   2^15
% m:7   max group order:8   group size:2097152 2^21

aux = eye(m);
for i = 0:2^(m*(m-1)/2)-1
    a = dec2bin(i) - '0';
    empty = m*(m-1)/2 - size(a,2);
    if empty
        a = [zeros(1,empty) a];
    end
    ind = [];
    for j = 2:m
        ind = [ind j+m*(j-2):m*(j-1)];
    end
    A = aux;
    A(ind) = a;
    A = A';
    for k = 0:m-2
        if mod(A^(2^k),2) == aux
            GROUP.order(i+1) = 2^k;
            GROUP.element(i+1).matrix = A;
            break
        elseif k == m-2
            GROUP.order(i+1) = 2^(m-1);
            GROUP.element(i+1).matrix = A;
        end
    end
    disp(['Done: ' num2str(i/(2^(m*(m-1)/2)-1))])
end

disp(['m:' num2str(m) '   max group order:' num2str(max(GROUP.order)) '   group size:' num2str(2^(m*(m-1)/2))])
    
    