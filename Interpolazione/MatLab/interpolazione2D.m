%zi= interp2(x, y, z, xi, yi, 'method'), dove il metodo può essere : linear, nearest, spline
%x e y sono i punti di interpolazione, la z è la funzione valutata nei
%punti. xi e yi sono i punti in cui voglio conoscere l interpolante


%Consideriamo un paraboloide di rotazione z=x^2 y y^2

%[XX,YY]=meshgrid(linspace(-2,2,50),linspace(-2,2,50)); %li dispone a griglia
%Z=XX.^2+YY.^2;
%zi=interp2(XX,YY,Z,-1.88,1.88,'linear'); %abbiamo visualizzato prima il coso
%figure();
%surf(Z);
%zi_s=interp2(XX,YY,Z,-1.88,1.88,'spline');
%zi_n=interp2(XX,YY,Z,-1.88,1.88,'nearest');

%IMMAGINI

A = imread('ngc6543a.jpg');
imshow(A)
imshow(A(1:570,10:600,:),'InitialMagnification','fit')
zoom(10)
title('Original Image, 10x Zoom') %facendo lo zoom si nota che abbiamo perso dei pixel. Per risolvere si può, esempio, usare l'interpolazione
F = griddedInterpolant(double(A)); %griglia dei pixel (i canali RGB vanno convertiti in double)
[sx,sy,sz] = size(A);
xq = (1:0.1:sx)';      % questa istruzione crea un vettore di valori da 1 a sx con passo 0.1
yq = (1:0.1:sy)';
zq = (1:sz)';        

figure()
F.Method = 'linear';
vq = uint8(F({xq,yq,zq}));       %rifaccio il cast per poter visualizzare
imshow(vq(1:5700,150:5900,:),'InitialMagnification','fit')
zoom(10)
title('Linear method')

figure()
F.Method = 'cubic';
vq = uint8(F({xq,yq,zq}));
imshow(vq(1:5700,150:5900,:),'InitialMagnification','fit')
zoom(10)
title('Cubic method')


figure()
F.Method = 'spline';
vq = uint8(F({xq,yq,zq}));
imshow(vq(1:5700,150:5900,:),'InitialMagnification','fit')
zoom(10)
title('Spline method')


