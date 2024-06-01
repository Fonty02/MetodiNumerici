A = sparse([4, 0, 0; 0, 4, 0; 0, 0, 4]);
disp("A")
disp(full(A));
b = [1; 2; 3];
x0 = rand(3, 1);

disp("A COME MATRICE SPARSA")
tic;
soluzione = Jacobi(A, b, x0);
tempo=toc;
disp('Jacobi')
disp(soluzione)
disp('tempo')
disp(tempo)

tic;
soluzione = SORForward(A, b, x0);
tempo=toc;
disp('SORForward')
disp(soluzione)
disp('tempo')
disp(tempo)


tic;
soluzione = SORBackward(A, b, x0);
tempo=toc;
disp('SORBackward')
disp(soluzione)
disp('tempo')
disp(tempo)


tic;
soluzione = SORSymmetric(A, b, x0);
tempo=toc;
disp('SORSymmetric')
disp(soluzione)
disp('tempo')
disp(tempo)


tic;
soluzione = A\b;
tempo=toc;
disp('Solver classico')
disp(soluzione)
disp('tempo')
disp(tempo)

A = full(A);

tic;
soluzione = Jacobi(A, b, x0);
tempo=toc;
disp('Jacobi')
disp(soluzione)
disp('tempo')
disp(tempo)

tic;
soluzione = SORForward(A, b, x0);
tempo=toc;
disp('SORForward')
disp(soluzione)
disp('tempo')
disp(tempo)


tic;
soluzione = SORBackward(A, b, x0);
tempo=toc;
disp('SORBackward')
disp(soluzione)
disp('tempo')
disp(tempo)


tic;
soluzione = SORSymmetric(A, b, x0);
tempo=toc;
disp('SORSymmetric')
disp(soluzione)
disp('tempo')
disp(tempo)


tic;
soluzione = A\b;
tempo=toc;
disp('Solver classico')
disp(soluzione)
disp('tempo')
disp(tempo)


disp("MATRICE RANDOMICA A SPARSA GRANDE 10000x10000 con 90% di zeri")
A=sprand(10000,10000,0.9);
b=rand(10000,1);
x0=rand(10000,1);

%{
Il metodo di Jacobi potrebbe dare errori in caso di matrici singolari
tic;
soluzione1 = Jacobi(A, b, x0);
tempo=toc;
disp('Jacobi')
disp('tempo')
disp(tempo)
Jacobi
%}
tic;
soluzione = A\b;
tempo=toc;
disp('Solver classico')
disp('tempo')
disp(tempo)

%disp('Soluzioni uguali?')
%disp(isequal(soluzione, soluzione1))

disp("STESSA MATRICE A MA PIENA")
A=full(A);
%{

tic;
soluzione1 = Jacobi(A, b, x0);
tempo=toc;
disp('Jacobi')
disp('tempo')
disp(tempo)
%}
tic;
soluzione = A\b;
tempo=toc;
disp('Solver classico')
disp('tempo')
disp(tempo)

%disp('Soluzioni uguali?')
%disp(isequal(soluzione, soluzione1))

disp("MATRICE RANDOMICA A SPARSA GRANDE 10000x10000 con 40% di zeri")
A=sprand(10000,10000,0.4);
b=rand(10000,1);
x0=rand(10000,1);

%{
Il metodo di Jacobi potrebbe dare errori in caso di matrici singolari
tic;
soluzione1 = Jacobi(A, b, x0);
tempo=toc;
disp('Jacobi')
disp('tempo')
disp(tempo)
Jacobi
%}
tic;
soluzione = A\b;
tempo=toc;
disp('Solver classico')
disp('tempo')
disp(tempo)

%disp('Soluzioni uguali?')
%disp(isequal(soluzione, soluzione1))

disp("STESSA MATRICE A MA PIENA")
A=full(A);
%{

tic;
soluzione1 = Jacobi(A, b, x0);
tempo=toc;
disp('Jacobi')
disp('tempo')
disp(tempo)
%}
tic;
soluzione = A\b;
tempo=toc;
disp('Solver classico')
disp('tempo')
disp(tempo)

%disp('Soluzioni uguali?')
%disp(isequal(soluzione, soluzione1))




disp("METODO DELLE POTENZE PER MATRICE DI GOOGLE SPARSA")
n = 1000;
A = sprand(n, n, 0.15);
u0 = rand(n, 1);
t1 = tic;
[lam, u, n_it, err] = Met_PotenzeGoogle(u0, A);
t2 = toc(t1);
fprintf('tempo: %f\n', t2);
fprintf('autovalore: %f\n', lam);
fprintf('iterate: %d\n', n_it);
fprintf('errore: %f\n', err);
k=1;
[v, lambda] = eigs(A, k, 'largestabs');

% Visualizza l'autovalore più grande
disp('Autovalore più grande:');
disp(lambda);






function [D, E, F] = split(A)
    D = spdiags(diag(A), 0, size(A, 1), size(A, 2));
    E = tril(A, -1);
    F = triu(A, 1);
end

function x = Jacobi(A, b, x0, tol, max_iter)
    if nargin < 4
        tol = 1e-15;
    end
    if nargin < 5
        max_iter = 5000;
    end
    [D, E, F] = split(A);
    M = D;
    N = -(E + F);
    it = 0;
    stop = false;
    while it < max_iter && ~stop
        x1 = M \ (N * x0 + b);
        if norm(x1 - x0) < tol
            stop = true;
            break;
        end
        x0 = x1;
        it = it + 1;
    end
    if stop
        x = x1;
    else
        disp('Il metodo non converge');
        x = [];
    end
end

function x = SORForward(A, b, x0, tol, max_iter, omega)
    if nargin < 4
        tol = 1e-15;
    end
    if nargin < 5
        max_iter = 5000;
    end
    if nargin < 6
        omega = 1;
    end
    [D, E, F] = split(A);
    M = D - omega * E;
    N = -(omega * F + (1 - omega) * D);
    b = omega * b;
    it = 0;
    stop = false;
    while it < max_iter && ~stop
        x1 = M \ (N * x0 + b);
        if norm(x1 - x0) < tol
            stop = true;
            break;
        end
        x0 = x1;
        it = it + 1;
    end
    if stop
        x = x1;
    else
        disp('Il metodo non converge');
        x = [];
    end
end

function x = SORBackward(A, b, x0, tol, max_iter)
    if nargin < 4
        tol = 1e-15;
    end
    if nargin < 5
        max_iter = 5000;
    end
    [D, E, F] = split(A);
    M = D - F;
    N = -E;
    it = 0;
    stop = false;
    while it < max_iter && ~stop
        x1 = M \ (N * x0 + b);
        if norm(x1 - x0) < tol
            stop = true;
            break;
        end
        x0 = x1;
        it = it + 1;
    end
    if stop
        x = x1;
    else
        disp('Il metodo non converge');
        x = [];
    end
end

function x = SORSymmetric(A, b, x0, tol, max_iter, omega)
    if nargin < 4
        tol = 1e-15;
    end
    if nargin < 5
        max_iter = 5000;
    end
    if nargin < 6
        omega = 1;
    end
    [D, E, F] = split(A);
    ME = D - omega * E;
    NF = -(omega * F + (1 - omega) * D);
    MF = D - omega * F;
    NE = -(omega * E + (1 - omega) * D);
    b = omega * b;
    it = 0;
    stop = false;
    while it < max_iter && ~stop
        x1 = ME \ (NF * x0 + b);
        x2 = MF \ (NE * x1 + b);
        if norm(x1 - x0) < tol
            stop = true;
            break;
        end
        x0 = x2;
        it = it + 1;
    end
    if stop
        x = x1;
    else
        disp('Il metodo non converge');
        x = [];
    end
end

function [lam, u0, n_it, err, approx] = Met_PotenzeNorm(u0, A, tol, it_max)
    if nargin < 4
        tol = 1e-15;
    end
    if nargin < 5
        it_max = 5000;
    end
    n_it = 0;
    u1 = A * u0;
    u1 = u1 / norm(u1);
    lam0 = u1' * (A * u1) / (u1' * u1);
    u0 = u1;
    approx = lam0;
    err = 1;
    while err > tol && n_it < it_max
        u1 = A * u0;
        u1 = u1 / norm(u1);
        lam = u1' * (A * u1) / (u1' * u1);
        approx = [approx, lam];
        err = abs(lam - lam0) / (1 + abs(lam));
        lam0 = lam;
        u0 = u1;
        n_it = n_it + 1;
    end
end

function [lam, u0, n_it, err] = Met_PotenzeGoogle(u0, A, tol, it_max, alfa)
    if nargin < 4
        tol = 1e-15;
    end
    if nargin < 5
        it_max = 5000;
    end
    if nargin < 6
        alfa = 0.85;
    end
    n = size(A, 1);
    n_it = 0;
    u1 = alfa * A * u0 + ((1 - alfa) / n) * sum(u0);
    u1 = u1 / norm(u1);
    lam0 = u1' * (A * u1) / (u1' * u1);
    u0 = u1;
    err = 1;
    while err > tol && n_it < it_max
        u1 = alfa * A * u0 + ((1 - alfa) / n) * sum(u0);
        u1 = u1 / norm(u1);
        lam = u1' * (A * u1) / (u1' * u1);
        err = abs(lam - lam0) / (1 + abs(lam));
        lam0 = lam;
        u0 = u1;
        n_it = n_it + 1;
    end
end