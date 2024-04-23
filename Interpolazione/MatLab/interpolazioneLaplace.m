ref=imread('peppers.png');
ref=rgb2gray(ref); %rendo grigia l immagine
ref=double(ref);
[m,n]=size(ref);

%CREO UNA MASCHERA PER TOGLIERE ALCUNE INFORMAZIONI
P=rand(m,n);
P(P<=0.5)=1;
P(P~=1)=0; %diverso da 1
Rimossa=ref.*P;
[i,j]=find(Rimossa==0);
indici = [i, j];

for k = 1:length(indici)
    row=indici(k,1);
    col=indici(k,2);
    if row == 1 && col == 1 % angolo in alto a sinistra
        Rimossa(row, col) = (Rimossa(row + 1, col) + Rimossa(row, col + 1)) / 2;
    elseif row == 1 && col == n % angolo in alto a destra
        Rimossa(row, col) = (Rimossa(row + 1, col) + Rimossa(row, col - 1)) / 2;
    elseif row == m && col == 1 % angolo in basso a sinistra
        Rimossa(row, col) = (Rimossa(row - 1, col) + Rimossa(row, col + 1)) / 2;
    elseif row == m && col == n % angolo in basso a destra
        Rimossa(row, col) = (Rimossa(row - 1, col) + Rimossa(row, col - 1)) / 2;
    elseif row == 1 || row == m % prima riga o ultima riga
        Rimossa(row, col) = (Rimossa(row, col - 1) + Rimossa(row, col + 1)) / 2;
    elseif col == 1 || col == n % prima colonna o ultima colonna
        Rimossa(row, col) = (Rimossa(row - 1, col) + Rimossa(row + 1, col)) / 2;
    else % caso generale
        Rimossa(row, col) = (Rimossa(row - 1, col) + Rimossa(row + 1, col) + Rimossa(row, col - 1) + Rimossa(row, col + 1)) / 4;
    end
end

 figure()
 montage({uint8(ref),uint8(Rimossa)})
 title("Originale e LaplaceFattoBene")









%imshow(uint8(Rimossa))  %ho tolto in maniera sparsa il 25% dei pixel
 Mod=ref; %manteniamo il bordo
 Mod(2:end-1,2:end-1)=Rimossa(2:end-1,2:end-1);
% figure()
% montage({uint8(ref),uint8(Mod)})
% title("Originale e perturbata")

 %uso laplace male
 Recon=zeros(m,n);
 for i=2:m-2
     for j=2:n-2
         Recon(i,j)=(Mod(i-1,j)+Mod(i+1,j)+Mod(i,j-1)+Mod(i,j+1))/4;  %qui sto sostituendo tutti i pixel, non solo quelli mancanti
     end
 end

 Assembly=zeros(m,n);
 Assembly(1,:)=ref(1,:);
 Assembly(end,:)=ref(end,:);
 Assembly(:,1)=ref(:,1);
 Assembly(:,end)=ref(:,end);
 Assembly(2:end-1,2:end-1)=Recon(2:end-1,2:end-1);
 figure();
 montage({uint8(ref),uint8(Assembly)})
 title("Originale e Assembly")

% con le "tecniche gamma" posso "ingannare l occhio umano"
% Assembly=Assembly.^1.2;   %quando ci sono molti pixel neri vado a elevare per >1, per rendere piu scura <1
% figure();
% montage({uint8(ref),uint8(Assembly)})
% title("Originale e gamma")

