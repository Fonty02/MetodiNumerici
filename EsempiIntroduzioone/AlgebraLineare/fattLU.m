function [L, U] = factorizationLU(A)
    [m, n] = size(A);
    A = double(A); % non è necessario copiare la matrice in MATLAB
    tol = 1e-15;
    L = eye(n); % prendo la matrice identica di ordine n
    L=double(L);
    for k = 1:n-1
        for i = k+1:n
            mik = -A(i, k) / A(k, k);
            for j = k + 1:n
                A(i, j) = A(i, j) + mik * A(k, j);
            end
            L(i, k) = -mik;
        end
    end
    U = triu(A);
end