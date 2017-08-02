function m2 = fit_power_law(m)
    % regress best params
    [alpha,beta] = find_power_law(m);

    n = length(m);
    m2 = zeros(n, n);
    for i=1:n
        for j=1:n
            if i ~= j
                log_d = log(abs(i-j));
                m2(i,j) = exp(alpha*log_d + beta) - 1;
            end
        end
    end


end

function [alpha, beta] = find_power_law(m)
    n = length(m);
    n_entries = nchoosek(n, 2);
    A = ones(n_entries, 2);
    b = ones(n_entries);
    counter = 1;
    for i=1:n
        for j=i+1:n
            A(counter, 1) = log(abs(j-i));
            b(counter) = log(1 + m(i, j));
            counter = counter + 1;
        end
    end
    assert(counter == n_entries + 1);
    res = A\b;
    alpha = res(1);
    beta = res(2);
    
end
    
