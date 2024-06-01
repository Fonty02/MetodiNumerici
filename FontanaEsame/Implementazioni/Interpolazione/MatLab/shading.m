v1=[96,147,141];
v2=[126,207,242];
v3=[38,245,202];
v4=[194,18,89];

img=zeros(100,100,3);
img(1,1,:)=v1;
img(1,100,:)=v2;
img(100,1,:)=v3;
img(100,100,:)=v4;

%Interpolazione orizzontale
tic;
a = img(1, 1, :) .* ((1:100) -100)/(1-100) + img(100, 1, :) .* ((1:100) - 1)/(100-1);
b = img(1, 100, :) .* ((1:100) - 100)/(1-100) + img(100, 100, :) .* ((1:100) - 1)/(100-1);
%Interpolazione verticale -> fare (1:100)' perchè in questo modo si ottiene un vettore colonna
img = a .* (((1:100)' - 100) / (1 - 100)) + b .* (((1:100)' - 1) / (100 - 1));
tempo=toc;

imshow(uint8(img));
title(['Tempo di esecuzione: ', num2str(tempo), ' secondi']);
