%VISUALIZZARE LA SPARSITA'
%x=rand([1,10]);
%y=rand([1,10]);
%figure();
%voronoi(x,y)
%axis equal


%ESEMPIO MATLAB
x = -3 + 6*rand(50,1);
y = -3 + 6*rand(50,1);
v = sin(x).^4 .* cos(y);
F = scatteredInterpolant(x,y,v);
[xq,yq] = meshgrid(-3:0.1:3);

F.Method = 'nearest';          %usa voronoi -> infatti ci saranno zone "piatte", quelle con i valori costanti
vq1 = F(xq,yq);
plot3(x,y,v,'mo')
hold on
mesh(xq,yq,vq1)
title('Nearest Neighbor')
legend('Sample Points','Interpolated Surface','Location','NorthWest')

F.Method = 'linear';
vq2 = F(xq,yq);
figure
plot3(x,y,v,'mo')
hold on
mesh(xq,yq,vq2)
title('Linear')
legend('Sample Points','Interpolated Surface','Location','NorthWest')


F.Method = 'natural';
vq3 = F(xq,yq);
figure
plot3(x,y,v,'mo')
hold on
mesh(xq,yq,vq3)
title('Natural Neighbor')
legend('Sample Points','Interpolated Surface','Location','NorthWest')