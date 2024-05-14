A=imread('cameraman.tif');
%imshow(A);

sz=size(A);
xg=1:sz(1);
yg=1:sz(2);
F=griddedInterpolant({xg,yg},double(A));
h=0.1;

%Calcolo del gradiente
grad_1=(F({xg+h,yg})-2.*F({xg,yg})+F({xg,yg+h}))./h;
grad_2=(2.*F({xg,yg})-F({xg-h,yg})-F({xg,yg-h}))./h;
grad=[grad_1(1:end-1,:);grad_2(end,:)]; %tutte le righe tranne ultima in avanti, ultima indietro
grad=[grad(:,1:end-1) grad_2(:,end)]; %tutte le colonne tranne ultima in avanti, ultima indietro

%imshow(uint8(grad));

[counts,x]=imhist(uint8(grad));
T=otsuthresh(counts);
BW=imbinarize(grad,T);
%imshow(BW);


%CREO CON SVD
A=double(A);
[U,S,V]=svd(A);

n=50;
Rec_100=U(:,1:n)*S(1:n,1:n)*V(:,1:n)';
figure();
imshow(uint8(Rec_100));

Diff_i=abs(A-Rec_100);
F=griddedInterpolant({xg,yg},Diff_i);
grad_1=(F({xg+h,yg})-2.*F({xg,yg})+F({xg,yg+h}))./h;
grad_2=(2.*F({xg,yg})-F({xg-h,yg})-F({xg,yg-h}))./h;
grad=[grad_1(1:end-1,:);grad_2(end,:)]; %tutte le righe tranne ultima in avanti, ultima indietro
grad=[grad(:,1:end-1) grad_2(:,end)]; %tutte le colonne tranne ultima in avanti, ultima indietro
imshow(uint8(grad));
imshow(uint8(Diff_i));

Scaled=(Diff_i-min(min(Diff_i)))./(max(max(Diff_i))-min(min(Diff_i)));
