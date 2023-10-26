function discos_gersgorin(A)
%% discos_gersgorin(A) grafica los discos de Gersgorin por filas (azul) y columnas (rojo) de la matriz A.
%
%  Parámetros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%
%    Salida:
%           Gráfica de los discos de Gersgorin de A.
%
%
%  Compute los datos necesarios para iniciar la función.
%%
%%
B=abs(A); %Mariz de valor absoluto
n=size(A,1); %Número de filas de A

R=zeros(1,n);%Vector de radios por filas (Rows)
C=zeros(1,n);%Vector de radios por columnas (Columns)

X=zeros(1,n); %Centro X
Y=zeros(1,n); %Centro Y
for i=1:n;
     R(i)=sum(B(i,:))-B(i,i);
     C(i)=sum(B(:,i))-B(i,i);
     hold on
     X(i)=real(A(i,i));
     Y(i)=imag(A(i,i));
end;
circles(X,Y,R,'edgecolor','b','facecolor',[0.4258 0.4805 0.8320],...
    'facealpha',.4); %Grafica los discos por filas (rojo).
 hold on
circles(X,Y,C,'edgecolor','r','facecolor',[0.8594 0.2578 0.4180],...
    'facealpha',.4); %Grafica los discos por columnas (azul).

