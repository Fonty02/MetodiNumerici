x = [1, 2, 3, 4];
y = [2, 3, 5, 2];
xx = linspace(1, 4, 100);

tic; % Avvia il cronometro
yy = splineQuadratica(x, y, xx);
tempo = toc; % Arresta il cronometro

figure;
plot(x, y, 'ro', 'DisplayName', 'Punti dati');
hold on;
plot(xx, yy, 'b-', 'DisplayName', 'Interpolazione spline quadratica');
xlabel('x');
ylabel('y');
title(['Interpolazione spline quadratica con tempo di esecuzione: ', num2str(tempo)]);
legend;
grid on;
hold off;


function yy=splineQuadratica(x,y,xx)
    n=length(x);
    matrix=zeros(n+1,n+1);
    coeffienti=zeros(n+1,1);   %i coefficienti sono cosi organizzati-> b1 .. bn c2 .. cn
    termini_noti=zeros(n+1,1);
    b=zeros(n-1,1);
    c=zeros(n-1,1);


    %inserisco nella matrice i termini noti delle condizioni di continuità della funzione
    for i=1:n-1
        hi=x(i+1)-x(i);
        termini_noti(i)=y(i+1)-y(i);
        matrix(i,i)=hi;
        matrix(i,i+2)=hi^2;
    end
    matrix(1,3)=0;
    %inserisco nella matrice i termini noti delle condizioni di continuità della derivata prima
    j=1;
    for i=n:n+1
        matrix(i,j)=1;
        matrix(i,j+1)=-1;
        matrix(i,j+2)=2*(x(j+1)-x(j));
        j=j+1;
    end
    matrix(n,3)=0;
    coeffienti=matrix\termini_noti;
    for i=1:n-1
        b(i)=coeffienti(i);
    end
    j=2;
    for i=n:n+1
        c(j)=coeffienti(i);
        j=j+1;
    end
    yy=zeros(1,length(xx));
    for i=1:n-1
        index=find(xx>=x(i) & xx<=x(i+1));
        yy(index)=y(i)+b(i)*(xx(index)-x(i))+c(i)*(xx(index)-x(i)).^2;
    end
end
    
    