function [Unitaria,Triangular_Superior]=schur2(A)
%% descomposicion_schur(A) encuentra la descomposición aproximada de Schur de una matriz A (Usando el método de las potencias).
%
%  Parámetros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%          
%    Salida:
%           Unitaria: Matriz unitaria compleja (nxn).
%           Triangular_Superior: Matriz triangular superior compleja (nxn).
%
%  De modo que U'*A*U=T.
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
%A=[-8 -3 -13 -3; 30 11 38 10 ; 0 0 3 0 ; -3 -1 -5 0];
n=size(A,1);
Unitaria=eye(n); %Se almacena la matriz Unitaria
A1=A;
for k=1:n
    lambda=max(eig(A1));
    W=Vector_propio(A1,lambda);
    V=eye(n+1-k); %Matriz a ortonormalizar.  
    V(:,1)=W; %Vector propio asociado al primer valor propio.
    [U]=proceso_gram_schmidt(V); % U matriz ortogonal.
    C=suma_directa(eye(k-1),U); % Aumentar U para que su tamaño sea nxn.
    Unitaria=Unitaria*C; %Factores de la matriz unitaria en la descomposición.
    B=U'*A1*U; %Matriz por bloques [lambda1]+[A1].
    A1=B(2:n+1-k,2:n+1-k); %Tomar el bloque [A1].
end
Triangular_Superior=conj(Unitaria)'*A*Unitaria;