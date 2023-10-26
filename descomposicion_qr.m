function [Q,R]=descomposicion_qr(V)
%% descomposicion_qr(V) encuentra la descomposición QR de una matriz A.
%Mo Chen (2020). Gram-Schmidt orthogonalization 
%(https://www.mathworks.com/matlabcentral/fileexchange/55881-gram-schmidt-orthogonalization), 
%MATLAB Central File Exchange. Retrieved November 17, 2020.
%
%  Parámetros: 
%
%    Entrada:
%           V: Matriz cuadrada compleja (nxn) (cuyas columnas son los
%           vectores a ortonormalizar).    
%
%    Salida:
%           Q: Matriz cuadrada ortogonal compleja (nxn). 
%           R: Matriz triangular superior compleja (nxn).
%
%  De modo que A=Q*R.
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
n=size(V,1);
R=zeros(n); %Se almacena la matriz triangular superior.
Q=zeros(n); %Se almacena la matriz ortogonal.
R(1,1)=norm(V(:,1));
Q(:,1)=V(:,1)/R(1,1); %Normalizar primer vector.
for k=2:n
R(1:k-1,k)=Q(:,1:k-1)'*V(:,k); 
Q(:,k)=V(:,k)-Q(:,1:k-1)*R(1:k-1,k); %Restar la proyección.
R(k,k)=norm(Q(:,k));
Q(:,k)=Q(:,k)/R(k,k); %Normalizar el k-ésimo vector.
end