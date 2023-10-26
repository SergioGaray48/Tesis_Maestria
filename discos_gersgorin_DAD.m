function discos_gersgorin_DAD(A,D)
%% discos_gersgorin_DAD(A,D) grafica los discos de Gersgorin por filas (azul) y columnas (rojo) de la matriz D^-1AD.
%
%  Par�metros:
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%           D: Matriz cuadrada diagonal compleja (nxn).
%    Salida:
%           Gr�fica de los discos de Gersgorin de D^-1AD.
%
%
%  Compute los datos necesarios para iniciar la funci�n.
%%
%%
A1=D^(-1)*A*D;
discos_gersgorin(A1)
