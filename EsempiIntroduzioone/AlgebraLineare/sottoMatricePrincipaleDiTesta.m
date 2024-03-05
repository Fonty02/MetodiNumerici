function sm=sottoMatricePrincipaleDiTesta(A,k)
    [m,n]=size(A);
    if k>m || k>n || k<1
        disp("ERRORE");
        sm=zeros(m,n);
        return;
    end
    sm=zeros(k,k);
    for i=1:k
        for j=1:k
            sm(i,j)=A(i,j);
        end
    end
end