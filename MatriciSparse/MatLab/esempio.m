A=sprandn(10,10,0.1) % crea matrice sparsa -> FORMATO COORDINATE      coordinate  valore
nnz(A)  %numero elementi NON nulli
B=full(A) %visualizzare come piena

%con tic, operazione, toc          prendere i tempi circa 10 volte e fare la media

%Per passare da piena a sparsa
Bs=sparse(B) %FORMATO COORDINATE      (devono essere degli zeri PRECISI)

spy(Bs)

whos %per vedere info