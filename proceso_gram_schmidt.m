function [U]=proceso_gram_schmidt(V)
%% proceso_gram_schmidt(V) ortonormaliza el conjunto de vectores columna de V.
%
%  Parámetros: 
%
%    Entrada:
%           V: Matriz cuadrada compleja (nxn) (cuyas columnas son los
%           vectores a ortonormalizar).    
%
%    Salida:
%           U: Matriz cuadrada compleja (nxn) cuyas columnas forman un
%           conjunto ortonormal.
%
%  De modo que U(:,i)'*U(:,j)=0 si i~=j y U(:,i)'*U(:,i)=1.
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
n = size(V,2); 
U(:,1) = V(:,1)/norm(V(:,1)); %Normalizar el primer vector.
for j = 2:n
    v = V(:,j); 
    u = v;
    for k = 1:j-1
        w = U(:,k);
        u = u - ((w'*v)/(w'*w))*w; %Restarle la proyección de v sobre w.
    end
    U(:,j) = u/norm(u); %Normalizar cada vector.
end

