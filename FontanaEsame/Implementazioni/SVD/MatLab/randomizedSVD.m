A = randn(10, 20);
l=5;
q=5;
[U,S,V] = randomizedSVD_v1(A,l);
A_ricostruita = U*S*V';
disp("NORMA DIFFERENZA = "+norm(A-A_ricostruita));
[U,S,V] = randomizedSVD_v2(A,l,q);
A_ricostruita = U*S*V';
disp("NORMA DIFFERENZA = "+norm(A-A_ricostruita));
[U,S,V] = randomizedSVD_v3(A,l,q);
A_ricostruita = U*S*V';
disp("NORMA DIFFERENZA = "+norm(A-A_ricostruita));



function [U, S, V] = randomizedSVD_v1(A, l)
    [m, n] = size(A);
    Omega = randn(n, l);
    Y = A * Omega;
    [Q, ~] = qr(Y);
    B = Q * A;
    [U_tilde, S, V] = svd(B);
    U = Q * U_tilde;
end



function [U,S,V] = randomizedSVD_v2(A,l,q)
    [m,n] = size(A);
    Omega = randn(n,l);
    Y0 = A*Omega;
    [Q0,~] = qr(Y0);
    for i = 1:q
        Yj = A' * Q0;
        [Qj,~] = qr(Yj);
        Y0 = A * Qj;
        [Q0,~] = qr(Y0);
    end
    B = Q0*A;
    [U_tilde,S,V] = svd(B);
    U = Q0*U_tilde;
end

function [U,S,V] = randomizedSVD_v3(A,l,q)
    [m,n] = size(A);
    Omega = randn(m,l);
    Y=A*A';
    for i = 1:q
        Y = Y*Y;
    end
    [Q, ~] = qr(Y);
    B = Q * A;
    [U_tilde, S, V] = svd(B);
    U = Q * U_tilde;
end



