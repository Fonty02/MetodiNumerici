 function [x1,x2,delta]=eq2(a,b,c)
 % calcola le radici di unequazione di secondo grado
 delta=b^2-4*a*c;
 x1=(-b-sqrt(delta))/(2*a);
 x2=(-b+sqrt(delta))/(2*a);

 