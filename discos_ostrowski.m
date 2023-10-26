function discos_ostrowski(A,alpha)
%% discos_ostrowski(A,alpha) grafica los discos de Ostrowski de la matriz A.
%
%  Parámetros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%           alpha: Valor en [0,1].
%    Salida:
%           Gráfica de los discos de Ostrowski de A.
%
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
B=abs(A); %Mariz de valor absoluto
T=size(A,1); %Número de filas de A

R=zeros(1,T);%Radio Filas (Rows)
C=zeros(1,T);%Radio Columnas (Columns)

X=zeros(1,T); %Centro X
Y=zeros(1,T); %Centro Y
for i=1:T;
     R(i)=sum(B(i,:))-B(i,i);
     C(i)=sum(B(:,i))-B(i,i);
     hold on
X(i)=real(A(i,i));
Y(i)=imag(A(i,i));
end;

r=(R.^(alpha)).*(C.^(1-(alpha))); %Radio de los discos de Ostrowski.

hold on

circles(X,Y,r,'edgecolor','g','facecolor',[0.2070 0.7461 0.5469],...
    'facealpha',.4); %Grafica los discos de Ostrowski