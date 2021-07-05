function fixed_data = fixnan(data)
    fixed_data = data ;
    [M, N] = size(data);
    for m=1:M
        for n=1:N
            if isnan(data(m,n))
                s = 0;
                c = 0;
                if m > 1
                    s = s + data(m-1,n);
                    c = c+1;
                end
                if m < M
                    s = s + data(m+1,n);
                    c = c+1;
                end
                if n > 1
                    s = s + data(m,n-1);
                    c = c+1;
                end
                if n < N
                    s = s + data(m,n+1);
                    c = c+1;
                end
                fixed_data(m,n) = s/c;
            end
        end
    end
end