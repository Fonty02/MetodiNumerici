x=[1 2 3 4];
f=[96 126 38 194];
a=interpolazioneNewton(x,f);
t=1:0.01:4; %valori di t per cui calcolare il polinomio
polynomial = pol(a,x);
y = polynomial(t);
plot(t,y);
hold on;
scatter(x,f);
title('Interpolazione di Newton',a);
xlabel('x');
ylabel('f(x)');
legend('Polinomio interpolante','Punti di interpolazione');
hold off;



function a=interpolazioneNewton(x,f)
    diff_div=f;
    a(1)=diff_div(1);
    for i=2:4
        diff_div=(diff_div(2:end)-diff_div(1:end-1))./(x(i:end)-x(1:end-(i-1))); 
        a(i)=diff_div(1);
    end
end