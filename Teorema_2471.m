function [A_eps]=Teorema_2471(A,epsilon)
%% Teorema_2471(A,eps) calcula una matriz A_eps diagonalizable de modo que norm(A-A_eps,fro)^2<eps.
%
%  Par�metros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%           epsilon: Determina que tan lejos se quiere la matriz diagonalizable
%           A_eps de A.
%           u: Vector inicial complejo (nx1) (para aplicar m�todo de las
%           potencias).
%           it_max: N�mero m�ximo de iteraciones (para aplicar m�todo de las
%           potencias).
%           tol: Tolerancia en el error (para aplicar m�todo de las
%           potencias).
%
%    Salida:
%           A_eps: Matriz diagonalizable (nxn)
%
%  De modo que norm(A-A_eps,fro)^2<eps.
%
%  Compute los datos necesarios para iniciar la funci�n.
%%
n=size(A,1); % Tama�o de A.
%[U,~]=descomposicion_schur_potencias(A,u,it_max,tol);
[U,~]=schur(A); % Se obtiene la matriz unitaria en la descomposici�n de Schur.
cota=sqrt(epsilon/n); % Cota para las entradas diagonales de E.
ei=cota*rand(1,n); % Entradas diagonales de E (aleatorias).
E=diag(ei); % Matriz diagonal que perturba a T.
A_eps=A+U'*E*U; % Matriz de salida.
norm(A_eps-A); % Comprobaci�n: norm(A_eps-A)<eps.