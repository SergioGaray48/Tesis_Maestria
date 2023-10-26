function [S_eps,T_eps]=Teorema_2472(A,epsilon)
%% Teorema_2472(A,eps) calcula una matriz T_eps triangular superior similar a A v�a S_eps, de modo que el m�dulo de las entradas de T_eps sobre la diagonal principal son menores a eps. 
%
%  Par�metros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%           epsilon: Determina que tan peque�as se quieren las entradas sobre
%           la diagonal principal de T_eps.
%           u: Vector inicial complejo (nx1) (para aplicar m�todo de las
%           potencias).
%           it_max: N�mero m�ximo de iteraciones (para aplicar m�todo de las
%           potencias).
%           tol: Tolerancia en el error (para aplicar m�todo de las
%           potencias).
%
%    Salida:
%           S_eps: Matriz invertible (nxn)
%           T_eps: Matriz triangular superior (nxn)
%
%  De modo que S_eps^-1*A*S_eps=T_eps y se tiene que abs(T(i,j))<eps para i<j.
%
%  Compute los datos necesarios para iniciar la funci�n.
%%
%%
n=size(A,1); % Tama�o de A.
%[U,T]=descomposicion_schur_potencias(A,epsilon,u,it_max,tol); % Se obtiene la descomposici�n de Schur (Potencias).
[U,T]=schur(A,'complex');
T1=triu(T,1); %Parte triangular superior de T sin diagonal
t=max(max(T1)); 
v=epsilon.^(0:n-1);
D_eps=diag(v);
if t<=1;
    S_eps=U*D_eps;
else 
    w=t.^-(0:n-1);
    D_t=diag(w);
    S_eps=U*D_t*D_eps;
end
T_eps=S_eps^-1*A*S_eps;