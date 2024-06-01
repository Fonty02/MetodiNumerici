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
 title("Originale e Laplace")
