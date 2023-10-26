function raiz =raiz_cercana(a,v)
%% RaizCuad(a,v) calcula la raiz del polinomio cuadr�tico determinado por v m�scercana al punto a.
%
%  Par�metros:
%
%    Entrada:
%           a: N�mero complejo 
%           v: Vector con los coeficientes de un polinomio cuadr�tico (1x3).
%           
%    Salida:
%           raiz: Raiz del polinomio determinado por v m�s cercana al punto a.
%
%
%  Compute los datos necesarios para iniciar la funci�n.
%%
discriminante = v(2)^2-4*v(1)*v(3);
    r1 = -2*v(3)/(v(2)+sqrt(discriminante));
    r2 = -2*v(3)/(v(2)-sqrt(discriminante));
   if  abs(a-r1)-abs(a-r2)>=0
     raiz = r2;
    else raiz = r1;
    end;