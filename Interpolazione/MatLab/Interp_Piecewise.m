% Interpolazione 1D mediante il comando interp1

t = [0 20 40 56 68 80 84 96 104 110];
v = [0 20 20 38 80 80 100 100 125 125];
tt = linspace(0,110);
vl = interp1(t,v,tt);
plot(t,v,'o',tt,vl,'LineWidth',1.5)


hold on
vn = interp1(t,v,tt,'nearest');
plot(tt,vn,'LineWidth',1.5)

vs = interp1(t,v,tt,'spline'); %Spline cubica naturale
plot(tt,vs,'LineWidth',1.5)

vh = interp1(t,v,tt,'pchip'); % Hermite cubica
plot(tt,vh,'LineWidth',1.5)

m = makima(t,v,tt); %Akima modificata
hold on
plot(tt,m,'LineWidth',1.5)

legend('Dati','lineare','nearest','spline','pchip','makima', Location='best')


% Costruzione spline quadratica 
dati = [3.0;4.5;7.0;9.0]; %nodi
f = [2.5,1.0,2.5,0.5]; % valori funzionali


%termini noti: 
tn = [f(2)-f(1); f(3)-f(2); f(4)-f(3); 0;0];

a = f;

%matrice del sistema lineare:

h = diff(dati); %vettore della spaziatura tra due nodi successivi

%matrice sistema lineare:
A =[h(1),0,0,0,0; 
      0, h(2), h(2)^2,0,0; 
      0,0,0,h(3), h(3)^2;
       1,-1,0,0,0;
       0,1,2*h(2),-1,0];

sol  = A\tn; % sol =  b_1,b_2, c_2,b_3,c_3

grado = 2;
coefficienti = zeros(length(h), grado+1);

coefficienti(1,:) = [  0, sol(1),a(1)]; %c_1 ==0
coefficienti(2,:) = [ sol(3), sol(2), a(2)];
coefficienti(3,:) = [ sol(5),sol(4),a(3)];

% Creiamo una struttura di tipo "polinomio a tratti"
% pp = mkpp(breaks,coefs)

pp = mkpp(dati, coefficienti);

%valutiamo il polinomio creato:
xq = linspace(dati(1), dati(end)); %creo 100 punti di valutazione
yq = ppval(pp,xq);

figure()
plot(dati, f, 'o',xq, yq,'-b','LineWidth',1.5)
