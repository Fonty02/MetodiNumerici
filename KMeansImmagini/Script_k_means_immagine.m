
%Posizionare la cartella di lavoro dove si è scaricato l'immagine.

% Task: leggere l'immagine e identificare eventuali aree verdi

%Leggere l'immagine e memorizzarla in una variabile A
A = imread('pre_04_cambiamento_ro_cia_50cm.tif');
I = A;
figure()
imshow(A);
% come si può vedere dalla "workspace" A è un tensore di dimensione 2530 x 2990 x 3 

%Trasformare l'immagine in "double"
A = double(A);

%Prendere la banda del verde, cioè la seconda:
Verde = A(:,:,2);

%Riscalare i valori da [0, 255] -- > [0, 1]: 
% Usare una retta del tipo y - y0 = m(x-x0);
MM = max(max(Verde));
mm = min(min(Verde));
m = 1/(MM-mm);

New_verde = m*(Verde-MM)+1;

%Usare K-means per creare due gruppi
idx = kmeans(reshape(New_verde,2530*2990,1),2);

% "Colorare" i pixel di A secondo le etichette date come output da Kmeans:
S= labeloverlay(I,reshape(idx,2530,2990));
figure()
imshow(S)
