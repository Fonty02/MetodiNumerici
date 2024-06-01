dati = readtable('winequality-red.csv', 'Delimiter', ';');

rng(5); %Random seed per riproducibilità
cv = cvpartition(size(dati, 1), 'Holdout', 0.2); %Split 80-20% tra training e test set

% Training set
X_train = dati(cv.training, 1:11); % Prendi le colonne da 1 a 11 per le features di input
y_train = dati(cv.training, 12); % Prendi la colonna 12 per la feature di output

% Testing set
X_test = dati(cv.test, 1:11); 
y_test = dati(cv.test, 12); 

nd = size(X_train, 1); % numero dati di training

% Costruzione della matrice dei predittori per i dati di addestramento
A = ones(nd, 12); % Creazione di una matrice di tutti 1 con 12 colonne
%nelle coolonne da 2 a 12 inseriscoo i dati di X_train
A(:, 2:12) = table2array(X_train);
yest = table2array(y_train); % Dati noti

%Risoluzione del problema ai minimi quadrati (calcolo la QR economica)
tic;
[Q,R]=qr(A,'econ');

soluzione = R \ (Q' * yest);
tempo = toc;
fprintf('I coefficienti della regressione con QR sono:\n');
disp(soluzione);
fprintf('Il tempo impiegato per la risoluzione con QR è: %f\n', tempo);


%risoluzione del problema ai minimi quadrati con la SVD
tic;
[U, S, V] = svd(A, 0);
soluzione = V * (S \ (U' * yest));
tempo = toc;
fprintf('I coefficienti della regressione con SVD sono:\n');
disp(soluzione);
fprintf('Il tempo impiegato per la risoluzione con SVD è: %f\n', tempo);


%risoluzione con equazioni normali
tic;
soluzione = (A' * A) \ (A' * yest);
tempo = toc;
fprintf('I coefficienti della regressione con equazioni normali sono:\n');
disp(soluzione);
fprintf('Il tempo impiegato per la risoluzione con equazioni normali è: %f\n', tempo);


%TESTING -> predico i valori di y per i dati di test, li confronto con i valori noti e calcolo l'errore RMSE (Root Mean Squared Error)
nd = size(X_test, 1); % numero dati di test
A = ones(nd, 12); % Creazione di una matrice di tutti 1 con 12 colonne
A(:, 2:12) = table2array(X_test); % inserisco i dati di X_test

yest = A * soluzione; % predizione dei valori di y
y_test = table2array(y_test); % valori noti

errore = sqrt(mean((yest - y_test).^2)); % calcolo dell'errore RMSE
fprintf("L'errore RMSE è: %f\n", errore);
