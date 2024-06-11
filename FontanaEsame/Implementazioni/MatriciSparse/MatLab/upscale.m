
%IMMAGINI

A = imread('C:\Users\fonta\Desktop\sparkingzerosfondo.jpg');
F = griddedInterpolant(double(A)); %griglia dei pixel (i canali RGB vanno convertiti in double)
[sx,sy,sz] = size(A);
%resize to 1920x1080
xq = linspace(1,sx,1080);
yq = linspace(1,sy,1920);
zq = 1:3;        

size(A)
figure()
F.Method = 'spline';
vq = uint8(F({xq,yq,zq}));
imshow(vq,'InitialMagnification','fit')
size(vq)
%save the image
imwrite(vq,'C:\Users\fonta\Desktop\sparkingzerosfondo1920x1080.jpg');


