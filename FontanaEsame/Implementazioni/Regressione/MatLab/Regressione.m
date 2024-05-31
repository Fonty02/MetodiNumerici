dati = readtable('winequality-red.csv', 'Delimiter', ';');

rng(5); % Set random seed for reproducibility
cv = cvpartition(size(dati, 1), 'Holdout', 0.2); % Hold out 20% of the data for testing

% Training set
X_train = dati(cv.training, 1:11); % Select columns 1 through 11 for features
y_train = dati(cv.training, 12); % Select column 12 for labels

% Testing set
X_test = dati(cv.test, 1:11); % Select columns 1 through 11 for features
y_test = dati(cv.test, 12); % Select column 12 for labels

nd = size(X_train, 1); % numero dati di training

% Costruzione della matrice dei predittori per i dati di addestramento
A = ones(nd, 12); % Creazione di una matrice di tutti 1 con 12 colonne
%nelle coolonne da 2 a 12 inseriscoo i dati di X_train
A(:, 2:12) = table2array(X_train);
yest = table2array(y_train); % Dati noti

%Risoluzione del problema ai minimi quadrati (calcolo la QR economica)
[Q,R]=qr(A,'econ');

soluzione = R \ (Q' * yest);
fprintf('I coefficienti della regressione con QR sono:\n');
disp(soluzione);

%risoluzione del problema ai minimi quadrati con la SVD
[U, S, V] = svd(A, 0);
soluzione = V * (S \ (U' * yest));
fprintf('I coefficienti della regressione con SVD sono:\n');
disp(soluzione);



