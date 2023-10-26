function [C,poli_carac]=polinomio_caracteristico(A)
%% polinomio_caracteristico(A) computa el polinomio caracter�stico de la matriz A.
%
%  Par�metros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%           
%    Salida:
%           poli_carac: Polinomio caracteristico de A (grado n).
%
%
%  Compute los datos necesarios para iniciar la funci�n.
%%
%%
n=size(A,1);
x=zeros(1,n+1); % Vector con coordenadas de x
y=zeros(1,n+1); % Vector con coordenadas de y
for k=1:n+1
    x(k)=k-1; % Enteros del 0 a n.
    y(k)=det(x(k)*eye(n)-A); % Im�genes de x bajo el polinomio caracter�stico
end
[C,poli_carac]=interpolacion_lagrange(x,y);
