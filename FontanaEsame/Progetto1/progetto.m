%Preparazione di tutti i dati necessari
format long g
im1 = imread('Bern1.bmp');
im2 = imread('Bern2.bmp');
if (size(im1,3)>1)
    im1 = rgb2gray(im1);
end
if(size(im2,3)>1)
    im2 = rgb2gray(im2);
end
im1=double(im1);
im2=double(im2);
Diff_Mat=abs(im1-im2);
rank_diff=rank(Diff_Mat);
fprintf("RANK DELLA MATRICE DIFFERENZA=%d\n\n\n",rank_diff);
[U,S,V]=svd(Diff_Mat);
valori_singolari=diag(S);
fj = (valori_singolari(1:rank_diff).^2) / sum(valori_singolari(1:rank_diff).^2);

%PRIMO CRITERIO DI TRONCAMENTO
fprintf("SVD TRONCATA -> CRITERIO GUTTMAN-KEISER\n");
tic;
valori_singolari_troncati=valori_singolari(valori_singolari>1);
numero_valori_singolari_troncati=length(valori_singolari_troncati);
fprintf("NUMERO DI VALORI SINGOLARI MAGGIORI DI 1= %d\n",numero_valori_singolari_troncati);
U_troncato=U(:,1:numero_valori_singolari_troncati);
V_troncato=V(:,1:numero_valori_singolari_troncati);
S_troncato=S(1:numero_valori_singolari_troncati,1:numero_valori_singolari_troncati);
Diff_Mat_troncata=U_troncato*S_troncato*V_troncato';
tempo=toc;
fprintf("TEMPO IMPIEGATO PER IL TRONCAMENTO=%f\n",tempo);
[righe,colonne]=size(Diff_Mat_troncata);
fprintf("ERRORE DI TRONCAMENTO IN NORMA 2=%f\n",norm(Diff_Mat-Diff_Mat_troncata,2))
fprintf("ERRORE DI TRONCAMENTO IN NORMA FROBENIUS=%f\n\n\n",norm(Diff_Mat-Diff_Mat_troncata,'fro'));
Diff_Mat_Troncata_vettore=Diff_Mat_troncata(:);
X_min = min(Diff_Mat_Troncata_vettore);
X_max = max(Diff_Mat_Troncata_vettore);
Diff_Mat_Troncata_vettore = (Diff_Mat_Troncata_vettore - X_min) / (X_max - X_min);
soglie=[0.25,0.5,0.75,graythresh(Diff_Mat_Troncata_vettore)];
num_soglie = length(soglie);
figure();
for i = 1:num_soglie
    immagine_binaria = reshape(Diff_Mat_Troncata_vettore > soglie(i),righe,colonne);
    subplot(1, num_soglie+2, i);
    imshow(immagine_binaria);
    title(['Soglia ', num2str(soglie(i))]);
end
subplot(1, num_soglie+2, num_soglie+1);
imshow(Diff_Mat);
title('Immagine differenza');
subplot(1, num_soglie+2, num_soglie+2);
imshow(Diff_Mat_troncata);
title('Immagine differenza troncata');

%SECONDO CRIERIO DI TRONCAMENTO
fprintf("SVD TRONCATA -> CRITERIO DI TRONCAMENTO SCREENPLOT DI COTTEL\n");
tic;
figure();
plot(valori_singolari);
title('Screenplot');
xlabel('Indice');
ylabel('Valore singolare');
idx = findchangepts(valori_singolari,'MaxNumChanges', 4);
hold on;
plot(idx, valori_singolari(idx), 'ro', 'MarkerSize', 10);
legend('Valore singolare', 'Punto di gomito');
hold off;
numero_valori_singolari_troncati=idx(2);
fprintf("OSSERVANDO LO SCREENPLOT, DECIDO DI PRENDERE IL SECONDO GOMITO %d\n ",numero_valori_singolari_troncati);
U_troncato=U(:,1:numero_valori_singolari_troncati);
V_troncato=V(:,1:numero_valori_singolari_troncati);
S_troncato=S(1:numero_valori_singolari_troncati,1:numero_valori_singolari_troncati);
Diff_Mat_troncata=U_troncato*S_troncato*V_troncato';
tempo=toc;
fprintf("TEMPO IMPIEGATO PER IL TRONCAMENTO=%f\n",tempo);
[righe,colonne]=size(Diff_Mat_troncata);
fprintf("ERRORE DI TRONCAMENTO IN NORMA 2=%f\n ",norm(Diff_Mat-Diff_Mat_troncata,2));
fprintf("ERRORE DI TRONCAMENTO IN NORMA FROBENIUS=%f\n ",norm(Diff_Mat-Diff_Mat_troncata,'fro'));
fprintf("MINIMO NORMA 2 = %f\n",valori_singolari(numero_valori_singolari_troncati+1));
fprintf("MINIMO NORMA FROBENIUS = %f\n\n\n ",sqrt(sum(valori_singolari(numero_valori_singolari_troncati+1:rank_diff).^2)));
Diff_Mat_Troncata_vettore=Diff_Mat_troncata(:);
X_min = min(Diff_Mat_Troncata_vettore);
X_max = max(Diff_Mat_Troncata_vettore);
Diff_Mat_Troncata_vettore = (Diff_Mat_Troncata_vettore - X_min) / (X_max - X_min);
soglie=[0.25,0.5,0.75,graythresh(Diff_Mat_Troncata_vettore)];
num_soglie = length(soglie);
figure();
for i = 1:num_soglie
    immagine_binaria = reshape(Diff_Mat_Troncata_vettore > soglie(i),righe,colonne);
    subplot(1, num_soglie+2, i);
    imshow(immagine_binaria);
    title(['Soglia ', num2str(soglie(i))]);
end
subplot(1, num_soglie+2, num_soglie+1);
imshow(Diff_Mat);
title('Immagine differenza');
subplot(1, num_soglie+2, num_soglie+2);
imshow(Diff_Mat_troncata);
title('Immagine differenza troncata');



%TERZO CRITERIO DI TRONCAMENTO
fprintf("SVD TRONCATA -> CRITERIO DI TRONCAMENTO BASATO SULL'ENTROPIA\n");
tic;
entropiaMatrice = -sum(fj .* log(fj)) / log(rank_diff);
% Trova il primo indice dove la funzione cumulativa supera l'entropia
F = cumsum(fj);
b = find(F > entropiaMatrice, 1);
perc_ = 0.75;
numero_valori_singolari_troncati = fix(rank_diff * entropiaMatrice * perc_);
fprintf('Entropia = %f, suggerito k = %d\n', entropiaMatrice, numero_valori_singolari_troncati);
U_troncato=U(:,1:numero_valori_singolari_troncati);
V_troncato=V(:,1:numero_valori_singolari_troncati);
S_troncato=S(1:numero_valori_singolari_troncati,1:numero_valori_singolari_troncati);
Diff_Mat_troncata=U_troncato*S_troncato*V_troncato';
tempo=toc;
fprintf("TEMPO IMPIEGATO PER IL TRONCAMENTO=%f\n",tempo);
[righe,colonne]=size(Diff_Mat_troncata);
fprintf("ERRORE DI TRONCAMENTO IN NORMA 2 =%f\n ",norm(Diff_Mat-Diff_Mat_troncata,2));
fprintf("ERRORE DI TRONCAMENTO IN NORMA FROBENIUS =%f\n ",norm(Diff_Mat-Diff_Mat_troncata,'fro'));
fprintf("MINIMO NORMA 2 =%f\n",valori_singolari(numero_valori_singolari_troncati+1));
fprintf("MINIMO NORMA FROBENIUS =%f\n\n\n ",sqrt(sum(valori_singolari(numero_valori_singolari_troncati+1:rank_diff).^2)));
Diff_Mat_Troncata_vettore=Diff_Mat_troncata(:);
X_min = min(Diff_Mat_Troncata_vettore);
X_max = max(Diff_Mat_Troncata_vettore);
Diff_Mat_Troncata_vettore = (Diff_Mat_Troncata_vettore - X_min) / (X_max - X_min);
soglie=[0.25,0.5,0.75,graythresh(Diff_Mat_Troncata_vettore)];
num_soglie = length(soglie);
figure();
for i = 1:num_soglie
    immagine_binaria = reshape(Diff_Mat_Troncata_vettore > soglie(i),righe,colonne);
    subplot(1, num_soglie+2, i);
    imshow(immagine_binaria);
    title(['Soglia ', num2str(soglie(i))]);
end
subplot(1, num_soglie+2, num_soglie+1);
imshow(Diff_Mat);
title('Immagine differenza');
subplot(1, num_soglie+2, num_soglie+2);
imshow(Diff_Mat_troncata);
title('Immagine differenza troncata');



%QUARTO CRITERIO DI TRONCAMENTO
fprintf("SVD TRONCATA -> CRITERIO DI ENERGIA\n");
tic;
energia = sum(valori_singolari.^2);
perc_ = 0.9;
numero_valori_singolari_troncati = find(cumsum(valori_singolari.^2) / energia >= perc_, 1);
fprintf('Energia = %f, suggerito k = %d\n', energia, numero_valori_singolari_troncati);
fprintf("%d%% dell'energia totale=%f \n", perc_*100, perc_*energia);
fprintf('Energia troncata = %f\n', sum(valori_singolari(1:numero_valori_singolari_troncati).^2));
U_troncato=U(:,1:numero_valori_singolari_troncati);
V_troncato=V(:,1:numero_valori_singolari_troncati);
S_troncato=S(1:numero_valori_singolari_troncati,1:numero_valori_singolari_troncati);
Diff_Mat_troncata=U_troncato*S_troncato*V_troncato';
tempo=toc;
fprintf("TEMPO IMPIEGATO PER IL TRONCAMENTO=%f\n",tempo);
[righe,colonne]=size(Diff_Mat_troncata);
fprintf("ERRORE DI TRONCAMENTO IN NORMA 2 =%f\n",norm(Diff_Mat-Diff_Mat_troncata,2));
fprintf("ERRORE DI TRONCAMENTO IN NORMA FROBENIUS =%f\n ",norm(Diff_Mat-Diff_Mat_troncata,'fro'));
fprintf("MINIMO NORMA 2 =%f \n",valori_singolari(numero_valori_singolari_troncati+1));
fprintf("MINIMO NORMA FROBENIUS =%f\n\n\n",sqrt(sum(valori_singolari(numero_valori_singolari_troncati+1:rank_diff).^2)));
Diff_Mat_Troncata_vettore=Diff_Mat_troncata(:);
X_min = min(Diff_Mat_Troncata_vettore);
X_max = max(Diff_Mat_Troncata_vettore);
Diff_Mat_Troncata_vettore = (Diff_Mat_Troncata_vettore - X_min) / (X_max - X_min);
soglie=[0.25,0.5,0.75,graythresh(Diff_Mat_Troncata_vettore)];
num_soglie = length(soglie);
figure();
for i = 1:num_soglie
    immagine_binaria = reshape(Diff_Mat_Troncata_vettore > soglie(i),righe,colonne);
    subplot(1, num_soglie+2, i);
    imshow(immagine_binaria);
    title(['Soglia ', num2str(soglie(i))]);
end
subplot(1, num_soglie+2, num_soglie+1);
imshow(Diff_Mat);
title('Immagine differenza');
subplot(1, num_soglie+2, num_soglie+2);
imshow(Diff_Mat_troncata);
title('Immagine differenza troncata');


%QUINTO CRITERIO DI TRONCAMENTO
fprintf("SVD TRONCATA -> CRITERIO DI TRONCAMENTO K-MEANS ISOLATION FOREST\n");
tic;
log_fj = log(fj);
[idx,C] = kmeans(log_fj, 2);
cluster_informazione = find(idx == idx(1));
[forest,anomaly]=iforest(cluster_informazione);
if length(unique(anomaly))==1
    numero_valori_singolari_troncati = length(anomaly);
else
    numero_valori_singolari_troncati = find(diff(anomaly)==-1,1);
    %eseguo le differenze tra elementi successivi. Se trovo -1 significa che ho 1,0 e dunque un salto
end
fprintf("Numero di valori singolari troncati=%d\n",numero_valori_singolari_troncati);
U_troncato=U(:,1:numero_valori_singolari_troncati);
V_troncato=V(:,1:numero_valori_singolari_troncati);
S_troncato=S(1:numero_valori_singolari_troncati,1:numero_valori_singolari_troncati);
Diff_Mat_troncata=U_troncato*S_troncato*V_troncato';
tempo=toc;
fprintf("TEMPO IMPIEGATO PER IL TRONCAMENTO=%f\n",tempo);
[righe,colonne]=size(Diff_Mat_troncata);
fprintf("ERRORE DI TRONCAMENTO IN NORMA 2 =%f\n",norm(Diff_Mat-Diff_Mat_troncata,2));
fprintf("ERRORE DI TRONCAMENTO IN NORMA FROBENIUS =%f\n ",norm(Diff_Mat-Diff_Mat_troncata,'fro'));
fprintf("MINIMO NORMA 2 =%f \n",valori_singolari(numero_valori_singolari_troncati+1));
fprintf("MINIMO NORMA FROBENIUS =%f\n\n\n",sqrt(sum(valori_singolari(numero_valori_singolari_troncati+1:rank_diff).^2)));
Diff_Mat_Troncata_vettore=Diff_Mat_troncata(:);
X_min = min(Diff_Mat_Troncata_vettore);
X_max = max(Diff_Mat_Troncata_vettore);
Diff_Mat_Troncata_vettore = (Diff_Mat_Troncata_vettore - X_min) / (X_max - X_min);
soglie=[0.25,0.5,0.75,graythresh(Diff_Mat_Troncata_vettore)];
num_soglie = length(soglie);
figure();
for i = 1:num_soglie
    immagine_binaria = reshape(Diff_Mat_Troncata_vettore > soglie(i),righe,colonne);
    subplot(1, num_soglie+2, i);
    imshow(immagine_binaria);
    title(['Soglia ', num2str(soglie(i))]);
end
subplot(1, num_soglie+2, num_soglie+1);
imshow(Diff_Mat);
title('Immagine differenza');
subplot(1, num_soglie+2, num_soglie+2);
imshow(Diff_Mat_troncata);
title('Immagine differenza troncata');



