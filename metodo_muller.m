function [raiz,R,pasos,error_final]=metodo_muller(x,f,it_max,tol)
%% metodo_muller(x,f,it_max,tol) calcula una raiz compleja del polinomio que pasa por (x(i),f(i)) con un error especificado por tol.
%
%  Parámetros:
%
%    Entrada:
%           x: Vector con las coordenadas en x (1x3).
%           f: Vector con los coeficientes de un polinomio.
%           it_max: Número de iteraciones máximas.  
%           tol: Tolerancia del proceso (máximo error permitido).
%
%    Salida:
%           raiz: Raiz del polinomio f.
%           R: Almacena las iteraciones en el proceso de encontrar la raiz.
%           pasos: Pasos necesarios para alcanzar la raiz con un error
%           menor a la tolerancia.
%           error_final: Error final.
%
%  De modo que polyval(f,raiz)<tol.
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
   pasos=1;
   R=zeros(it_max,1); %Se almacenan las iteraciones en busca de la raíz.
   y=polyval(f,x);
  while pasos<=it_max && abs(y(3))> tol 
    [v,~]=interpolacion_muller(x,y);
    r=raiz_cercana(x(3),v);
    R(pasos)=x(3)+r;
    x(1)=x(2);
    x(2)=x(3);
    x(3)=x(3)+r;
    y=polyval(f,x);    
    pasos=pasos+1;
  end;
  error_final=y(3);
  raiz=R(pasos);