function raiz =raiz_cercana(a,v)
%% RaizCuad(a,v) calcula la raiz del polinomio cuadrático determinado por v máscercana al punto a.
%
%  Parámetros:
%
%    Entrada:
%           a: Número complejo 
%           v: Vector con los coeficientes de un polinomio cuadrático (1x3).
%           
%    Salida:
%           raiz: Raiz del polinomio determinado por v más cercana al punto a.
%
%
%  Compute los datos necesarios para iniciar la función.
%%
discriminante = v(2)^2-4*v(1)*v(3);
    r1 = -2*v(3)/(v(2)+sqrt(discriminante));
    r2 = -2*v(3)/(v(2)-sqrt(discriminante));
   if  abs(a-r1)-abs(a-r2)>=0
     raiz = r2;
    else raiz = r1;
    end;