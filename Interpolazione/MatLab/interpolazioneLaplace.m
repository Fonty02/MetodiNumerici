ref=imread('peppers.png');
ref=rgb2gray(ref); %rendo grigia l immagine
%imshow(ref)
ref=double(ref);
[m,n]=size(ref);

%CREO UNA MASCHERA PER TOGLIERE ALCUNE INFORMAZIONI
P=rand(m,n);
P(P<=0.75)=1;
P(P~=1)=0;
Rimossa=ref.*P;
%imshow(uint8(Rimossa))  %ho tolto in maniera sparsa il 25% dei pixel
Mod=ref; %manteniamo il bordo
Mod(2:end-1,2:end-1)=Rimossa(2:end-1,2:end-1);
figure()
montage({uint8(ref),uint8(Mod)})
title("Originale e perturbata")

%uso laplace
Recon=zeros(m,n);
for i=2:m-2
    for j=2:n-2
        Recon(i,j)=(Mod(i-1,j)+Mod(i+1,j)+Mod(i,j-1)+Mod(i,j+1))/4;  %qui sto sostituendo tutti i pixel, non solo quelli mancanti
    end
end
%bisognerebbe fare una lista dei pixel mancanti e agire solo su questi

Assembly=zeros(m,n);
Assembly(1,:)=ref(1,:);
Assembly(end,:)=ref(end,:);
Assembly(:,1)=ref(:,1);
Assembly(:,end)=ref(:,end);
Assembly(2:end-1,2:end-1)=Recon(2:end-1,2:end-1);
figure();
montage({uint8(ref),uint8(Assembly)})
title("Originale e Assembly")

%con le "tecniche gamma" posso "ingannare l occhio umano"
Assembly=Assembly.^1.2;   %quando ci sono molti pixel neri vado a elevare per >1, per rendere piu scura <1
figure();
montage({uint8(ref),uint8(Assembly)})
title("Originale e gamma")


%EX -> fai array con indici valori mancanti per correggere la Recon