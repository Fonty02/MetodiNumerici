%L applicazione è metterla in oggetti tassellati per fare altro

%Assegnamo dei colori in formato RGB ai 4 vertici. Aggiunta quarta colonna per la trasparenza
v1=[96,147,141];
v2=[126,207,242];
v3=[38,245,202];
v4=[194,18,89];

%Effettuiamo l'interpolazione bilineare tra i colori dei vertici per ottenere poi l'effetto di shading

%creazione di una matrice di 100x100x3 per contenere i colori
%img=zeros(100,100,3);
%img(1,1,:)=v1;
%img(1,100,:)=v2;
%img(100,1,:)=v3;
%img(100,100,:)=v4;
%for i=1:100
 %   for j=1:100
        %f(xi,y0)
        %a=img(1,1,:).*(i-100)/(1-100)+img(100,1,:).*(i-1)/(100-1);
        %f(xi,y1)
        %b=img(1,100,:).*(i-100)/(1-100)+img(100,100,:).*(i-1)/(100-1);
        %f(xi,yi)
        %img(i,j,:)=a.*(j-100)/(1-100)+b.*(j-1)/(100-1);
  %  end
%end
%imshow(uint8(img));



%Operazioni vettoriali (1 for)


img=zeros(100,100,3);
img(1,1,:)=v1;
img(1,100,:)=v2;
img(100,1,:)=v3;
img(100,100,:)=v4;
a = img(1, 1, :) .* ((1:100) -100)/(1-100) + img(100, 1, :) .* ((1:100) - 1)/(100-1);
b = img(1, 100, :) .* ((1:100) - 100)/(1-100) + img(100, 100, :) .* ((1:100) - 1)/(100-1);
for i=1:100
    img(i, :, :) = a .* (i-100)/(1-100) + b .* (i-1)/(100-1);
end
imshow(uint8(img));
