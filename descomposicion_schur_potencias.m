function [Unitaria,Triangular_Superior]=descomposicion_schur_potencias(A,u,it_max,tol)
%% descomposicion_schur(A) encuentra la descomposición aproximada de Schur de una matriz A (Usando el método de las potencias).
%
%  Parámetros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%           u: Vector inicial complejo (nx1) (para aplicar método de las
%           potencias). %Componentes no nulas.
%           it_max: Número máximo de iteraciones (para aplicar método de las
%           potencias).
%           tol: Tolerancia en el error (para aplicar método de las
%           potencias).
%
%    Salida:
%           Unitaria: Matriz unitaria compleja (nxn).
%           Triangular_Superior: Matriz triangular superior compleja (nxn).
%
%  De modo que U'*A*U=T (aproximadamente).
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
n=size(A,1);
Unitaria=eye(n); %Se almacena la matriz Unitaria
A1=A;
for k=1:n
    u=u(1:n+1-k);
    [~,W,~,~,~]=metodo_de_las_potencias(A1,u,it_max,tol); %Se encuentra vector propio asociado al valor propio dominante.
    V=eye(n+1-k); %Matriz a ortonormalizar.  
    V(:,1)=W; %Vector propio asociado al primer valor propio.
    [U]=proceso_gram_schmidt(V); % U matriz ortogonal.
    C=suma_directa(eye(k-1),U); % Aumentar U para que su tamaño sea nxn.
    Unitaria=Unitaria*C; %Factores de la matriz unitaria en la descomposición.
    B=U'*A1*U; %Matriz por bloques [lambda1]+[A1].
    A1=B(2:n+1-k,2:n+1-k); %Tomar el bloque [A1].
end
Triangular_Superior=Unitaria'*A*Unitaria;

