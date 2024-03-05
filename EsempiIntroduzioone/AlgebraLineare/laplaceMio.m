function res=laplaceMio(A)
[m,n]=size(A);
if n==1
    res=A(1,1);
else
    res=0;
    for j=1:n
        A1j=A(2:end,:);
        A1j=A1j(:,[1:j-1,j+1:end]);
        res=res+(-1)^(j-1) * A(1,j)*laplaceMio(A1j);
    end
end
end