function[val_prop,vec_prop,pasos,L,error_f]=metodo_de_las_potencias(A,u,it_max,tol)
%% metodo_de_las_potencias(A,u,tol) aplica el m�todo de las potencias a una matriz A, para encontrar el valor propio dominante y su respectivo vector propio normalizado.
%
%  Par�metros: 
%
%    Entrada:
%           A: Matriz cuadrada compleja (nxn).
%           u: Vector inicial complejo (nx1).
%           it_max: N�mero m�ximo de iteraciones.
%           tol: Tolerancia en el error.
%
%    Salida:
%           val_prop: Valor propio dominante de A
%           vec_prop: Vector propio unitario, asociado al valor propio dominante de A
%           pasos: N�mero de iteraciones para alcanzar la tolerancia del proceso           
%           L: Vector que almacena las iteraciones en el c�lculo del valor
%           propio dominante.
%           error_f: Error en el valor propio dominante.
%
%  De modo que A*(vec_prop)=(val_prop)*(vec_prop).
%
%  Compute los datos necesarios para iniciar la funci�n.
%%
%%
e=10; %e inicial para que e>tol.
%[p,~]=max(u); %Entrada de u con mayor valor absoluto.
%u=u/p; %Normalizar u (con la norma inf).
L=zeros(it_max,1); %Se almacenan las iteraciones para encontrar el valor propio dominante.
paso=1;
 while e>tol && paso<it_max 
   v=A*u; %Iteraciones m�todo de las potencias
   [p,~]=max(v); %Entrada de v con mayor valor absoluto.
   u=v/p;  %Normalizar v (con la norma inf).
   L(paso)=p; %Almacenar la sucesi�n que converge al valor propio dominante.
   e=norm(A*u-p*u); %Error entre iteraciones.
   paso=paso+1;
 end
pasos=paso-1;
L=L(1:pasos);
error_f=e; %Error final.
val_prop=L(pasos);
vec_prop=v;