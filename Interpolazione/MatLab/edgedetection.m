A = imread('cameraman.tif');
sz = size(A);
xg = 1:sz(1);
yg = 1:sz(2);
F = griddedInterpolant({xg,yg},double(A));

% Scelgo un incremento h:
 h = 0.1;

% Faccio sia le differenze in avanti che le differenze all'indietro 
% per l'approssimazione del gradiente 
grad_1 = (F({xg+h,yg}) - 2.*F({xg,yg}) + F({xg,yg+h}))./h;
grad_2 = (2.*F({xg,yg})-F({xg-h,yg})- F({xg,yg-h}))./h;

% Usiamo le differenze all'indietro per l'ultima riga e per l'ultima
% colonna:
grad = [grad_1(1:end-1,:);grad_2(end,:)];
grad = [grad(:,1:end-1),grad_2(:,end)];

figure()
imshow(uint8(grad))
title('Immagine gradiente')

%L'immagine gradiente ci sonsente di evidenziare i bordi.

%Per una migliore visualizzazione si usa un algoritmo di binarizzazione:
% Scelta una certa soglia (o più soglie) si genera un'immagine in bianco e
% nero. 

%Scaliamo l'immagine in [0-1]

Scaled = (grad - min(min(grad)))./(max(max(grad))-min(min(grad)));

Temp = Scaled;
Temp(Temp>=0.6) = 1;
Temp(Temp<0.6) = 0;
figure()
imshow(Temp)
title('Edge detection')



