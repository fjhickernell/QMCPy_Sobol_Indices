function y = t_value(C_j, C_d, m)
% T_VALUE gives the t-value between C_j and C_d when m is fixed
y = -1;
    for t = 0:m - 1
        if t == m - 1
            y = t;
            break
        end
        for r_j = 1:m - t
            r_d = m - t - r_j;
            M = [C_j(1:r_j,1:m); C_d(1:r_d,1:m)];
            M = mod(rref_toni(M),2);
            if ~any(M(end,:)) % Matrix M is not full rank (echelon matrix last row is 0)
                break
            elseif r_j == m - t
                y = t;
                break
            end
        end
        if y > -1
            break
        end
    end
end