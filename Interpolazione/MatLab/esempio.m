x=-3:3;
y=[-1 -1 -1 0 1 1 1];
xq1=-3:.01:3;
p=pchip(x,y,xq1); %HERMIT cubica per preservare la forma
s=spline(x,y,xq1);
m=makima(x,y,xq1);
figure();
plot(x,y,'o',xq1,p,'-',xq1,s,'-.',xq1,m,'--')
legend('Sample Points','pchip','spline','makima','Location','SouthEast')

%HERMIT E MAKIMA migliori quando c'è un po di linearità (vedi queste y)

pp=spline(x,y) %-> ritorna la struttura dove da: nodi (breaks), pezzi (intervalli), : ordine del polinomio (4 perchè spline 3), 
% coefficienti (nnhokpt)

%come valuto questo modello???

yq=ppval(pp,xq);

%interp1 -> funzione che ti fa scegliere la tipologia di interpolazione da
%usare

t=[]
v=[]
tt=[] %punti da interpolare
v1=interpo1(t,v,tt) %spline lienare
v2=interpo1(t,v,tt,'nearest') %un metodo basato su nearest neighbour
v3=interpo1(t,v,tt,'spline') %spline cubica
v4=interpo1(t,v,tt,'pchip') %hermit cubica
v5=interpo1(t,v,tt,'makima') %spline lienare
