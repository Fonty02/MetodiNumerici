function [A, U, fbar,M,ss1,ss2,N,RangeA] =svdDecompA


%OUTPUT:
% A matrice del nuovo gruppo di immagini 
% U matrice di vettori ortogonali che spannano R(A)
% fbar = immagine 'media'
% M = numero totale di pixels
% [ss1, ss2] = [righe,colonne] di ogni immagine
% N = numero totale di immagini presenti nel database
% RangeA = range di A

%percorso di locazione del trial set---aggiornare con la propria
%destinazione!!
myFolder = 'C:\Users\fonta\Desktop\MetodiNumerici\SVD\MatLab\faces94\faces94\female\9336923';

filePattern = fullfile(myFolder, '*.jpg');%cerco tutti e soli i files con estensione .jpg
theFiles = dir(filePattern);% lista i files nella current folder
N = length(theFiles);% numero totale di files presenti nella current folder


for i=1:N
    baseFileName = theFiles(i).name; %legge il nome del file in ordine lessicografico
    fullFileName = fullfile(myFolder, baseFileName);% costruisce il path al file contenuto in baseFileName
    fprintf('Now reading %s\n', fullFileName);
    
    % Trasformare l'immagine in '2D' passando alla scala di grigio per avere una MATRICE invece di un TENSORE
    fi = rgb2gray(imread(fullFileName));
    
    %Voglio visualizzare nella stessa finestra tutte le immagini presenti
    %nel database:
    subplot(ceil(sqrt(N)),ceil(sqrt(N)),i);
    figure(1); imshow(fi);
    
    %numero di righe e colonne di fi
    [ss1,ss2] = size(fi); %ss1=righe, ss2=colonne, 
    
    M = ss1*ss2; % DIMENSIONE

    S=zeros(M,N); %inizializzazione matrice delle facce
    
    fi = double(reshape(fi,M,1)); % ogni faccia diventa un vettore colonna
    
    S(:,i) = fi;%creazione matrice delle facce
  
end

fbar = mean(S,2); %immagine media

A = S-fbar; %matrice A -> ogni immagine sottrata la media

[U, ~, ~] = svd(A, 0); %siccome mi serve solo la matrice U per il range di A, non mi interessa calcolare SVD completa. 
RangeA = U(:,rank(A)); %Range di A
X = RangeA'*A;

%prendo una faccia della stessa persona (ovvero una colonna di S)
path_faccia_nota='C:\Users\fonta\Desktop\MetodiNumerici\SVD\MatLab\faces94\faces94\female\9336923\9336923.1.jpg';
faccia_nota = double(reshape(rgb2gray(imread(path_faccia_nota)),M,1));
%calcolo la proiezione della faccia nota sul range di A
proiezione = RangeA'*(faccia_nota-fbar);
%calcolo la distanza tra la proiezione e ogni colonna di X e prendo la più piccola
distanza_min = inf;
for i=1:N
    distanza = norm(proiezione-X(:,i));
    if distanza<distanza_min;
        distanza_min = distanza;
    end
end
figure()
imshow(reshape(uint8(faccia_nota+(fbar)),ss1,ss2)); %visualizzo la faccia nota
title(distanza_min)


%prendo una faccia di una persona DIVERSA
path_new_face = 'C:\Users\fonta\Desktop\MetodiNumerici\SVD\MatLab\faces94\faces94\female\anpage\anpage.3.jpg';
new_face = double(reshape(rgb2gray(imread(path_new_face)),M,1));
proiezione_new_face = RangeA'*(new_face-fbar);
distanza_min_new_face = inf;
for i=1:N
    distanza = norm(proiezione_new_face-X(:,i));
    if distanza<distanza_min_new_face;
        distanza_min_new_face = distanza;
    end
end
figure()
imshow(reshape(uint8(new_face+(fbar)),ss1,ss2)); %visualizzo la faccia nota
title(distanza_min_new_face)


end

