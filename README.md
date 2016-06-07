# -MOLER-COMPUTACI-N-NUM-RICA-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                                   %%%%%%%%%%%%%%
%%%%%%%%%%%%%             MOLER, COMPUTACIÓN NUMÉRICA           %%%%%%%%%%%%%%
%%%%%%%%%%%%%                                                   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                      %%%%%
%%%%%                       DEFINICIÓN DEL PROBLEMA:                       %%%%%
%%%%%   Implementar el método de cuadratura adaptativa usando la regla de  %%%%%     
%%%%%     Simpson como trabajado en clase. Complete el siguiente cuadro,   %%%%%
%%%%% donde los espacios en blanco corresponden al número de evaluaciones  %%%%%     
%%%%%     de la función necesarios para obtener la precisión deseada       %%%%%
%%%%%                                                                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                         VALORES DE ENTRADA:                     %%%%%
%%%%% fun: Es una función a la cual se desea calcular su integral     %%%%%     
%%%%% a,b: Son los extremos del intervalo donde deseamos integrar     %%%%%
%%%%%      la función fun                                             %%%%%   
%%%%% TOL: Es la tolerancia que deseamos tener en la aproximación     %%%%%   
%%%%%      De la integral                                             %%%%%   
%%%%%                                                                 %%%%%
%%%%%                         VALORES DE SALIDA:                      %%%%%
%%%%%                la función devuelve un vector [p N]              %%%%%     
%%%%% p: Es la aproximación a la integral de la función fun           %%%%%
%%%%%      la función fun                                             %%%%%   
%%%%% N: Es la cantidad de veces que se evaluo la función para        %%%%%   
%%%%%    tener el valor p                                             %%%%%   
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p N] = Cuadratura(fun,a,b,TOL)

APP = 0; 		                          %Aproximación para la integral
N = 3; 																%Cuantas iteraciones hacemos
i = 1;
TOL(i) = 10*TOL;
a(i) =  a;
h(i) = (b-a)/2;												%Divide el intervalo en dos
FA(i) = feval(fun,a); 								%Evalua la función en a
FC(i) = feval(fun,a + h(i)); 					%Evalua la función en a + h(i)
FB(i) = feval(fun,b);									%Evalua la función en a
S(i) = h(i)*(FA(i)+4*FC(i)+FB(i))/3;  %Aplica la Regla de Simpson
L(i) = 1;

while i>0 
	FD = feval(fun,a(i) + h(i)/2);       
	FE = feval(fun,a(i) + 3*h(i)/2);
	S1 = h(i)*(FA(i)+4*FD +FC(i))/6;
	
	S2 = h(i)*(FC(i)+4*FE +FB(i))/6;
	v1 = a(i);
	v2 = FA(i);
	v3 = FC(i);
	v4 = FB(i);
	v5 = h(i);
	v6 = TOL(i);
	v7 = S(i);
	v8 = L(i);
	
	i = i-1;
	if abs(S1+S2-v7)<v6
		APP = APP +(S1+S2);
	else
		i=i+1;                 %Los datos para la mitad derecha del intervalo
		a(i) = v1 + v5;
		FA(i) = v3;
		FC(i) = FE;
		FB(i) = v4;
		h(i) = v5/2;
		TOL(i) = v6/2;
		S(i) = S2;
		L(i) = v8+1;
		i = i+1;                %Los datos para la mitad izquierda del intervalo
		a(i) = v1;
		FA(i) = v2;
		FC(i) = FD;
		FB(i) = v3;
		h(i) = h(i-1);
		TOL(i) = TOL(i-1);
		S(i) = S1;
		L(i) = L(i-1);
	end
	N = N+2;
end
p = APP;                    %APP aproximaciones hechas dentro de TOL
