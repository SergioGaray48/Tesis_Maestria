function [C,polinomio_lagrange]=interpolacion_lagrange(x,y)
%% interpolacion_lagrange(x,y) calcula el polinomio que pasa por\ los valores especificados en (x,y).
%
%  Parámetros:
%
%    Entrada:
%           x: Vector con las coordenadas en x (1xn).
%           y: Vector con las coordenadas en f(x) (1xn).
%                
%    Salida:
%           C: Coeficientes del polinomio interpolador.
%           polinomio: Polinomio de grado n-1.
%
%  De modo que p(x(i))=y(i).
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
n=length(x); %Número de puntos a interpolar
L=zeros(n,n); % Se almacenan las bases polinómicas de Lagrange
for k=1:n
    V=1;
    for j=1:n
        if k~=j;
            V=conv(V,poly(x(j)))/(x(k)-x(j)); %poly(x(j))=[1 x(j)] %Se \ calculan las bases polinomicas de Lagrange.
        end
    end
    L(k,:)=V;
end
C=y*L; %Coeficientes del polinomio interpolador de Lagrange
polinomio_lagrange=poly2sym(C); %Polinomio interpolador
