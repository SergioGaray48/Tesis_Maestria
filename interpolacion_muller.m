function [coeficientes,polinomio]=interpolacion_muller(x,f)
%% interpolacion_muller(x,f) calcula el polinomio de grado 2 que pasa por (x(i),f(i)).
%
%  Parámetros:
%
%    Entrada:
%           x: Vector con las coordenadas en x (1x3).
%           f: Vector con las coordenadas en f(x) (1x3).
%           
%    Salida:
%           p: Polinomio de grado 2.
%
%  De modo que p(x(i))=f(i).
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
  coeficientes=[0;0;0];
 
    h1=x(2)-x(1);   
    h2=x(3)-x(2);
    d1=(f(2)-f(1))/h1;
    d2=(f(3)-f(2))/h2;
    
    coeficientes(1)= (d2-d1)/(x(3)-x(1));
    coeficientes(2)=d2+h2*coeficientes(1);
    coeficientes(3)=f(3);
%p(x)=v(1)(x-x(2))^2+v(2)(x-x(2))+v(3)    

polinomio=poly2sym(coeficientes(1)*conv(poly(x(2)),poly(x(2))))+poly2sym(coeficientes(2)*poly(x(2)))+poly2sym(coeficientes(3));
