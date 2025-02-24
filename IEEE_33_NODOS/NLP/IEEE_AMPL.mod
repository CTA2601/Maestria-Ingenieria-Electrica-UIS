#---------------------------------------------------------------------------
# Flujo de potencia no lineal formulado como un problema de
# optimización matemática
#---------------------------------------------------------------------------
# definir los conjuntos con el comando: set
#---------------------------------------------------------------------------
set Ob;						           # conjunto de barras
set Ol within Ob cross Ob;	           # conjunto de lineas del sistema
#---------------------------------------------------------------------------
# Parametros
#---------------------------------------------------------------------------
param R{Ol};		  # Resistencia de la línea
param X{Ol};		  # Reactancia de la línea
param Z2{Ol};		  # Z^2 = R^2 + X^2
param Imax{Ol};		  # Corriente máxima de la línea
param linea{Ol};	  # Número de la línea
param PD{Ob};		  # Potencia activa demandada
param QD{Ob};		  # Potencia reactiva demandada
param QC{Ob};		  # Potencia reactiva del condensador en la barra
param Ongd{Ob};		  # Indica si hay generación distribuida en la barra

param Vnom; 		         # Tensión nominal de la red
param Vmin;			         # Tensión mínima
param Vmax;			         # Tensión máxima
param Barra_SE;		         # Barra de la subestación
param Pmax_gd := 1;         # Potencia maxima que puede tener los GDs
param fp := 1;	             # Factor de potencia de la GD
param Nmax_gd := 10;         # Número máximo de GDs
param Ploss_max:= 0.202647;  # Límite máximo de pérdidas permitidas en la red
#---------------------------------------------------------------------------
# Se definen las variables con el comando var
#---------------------------------------------------------------------------
var I{Ol};			   # Corriente en las líneas
var P{Ol};		       # Potencia activa en las líneas
var Q{Ol};		       # Potencia reactiva en las líneas
var V{Ob};			   # Tensión en las barras
var Ps{Ob} >= -5;	   # Potencia activa generada en la subestación
var Qs{Ob} >= 0;	   # Potencia reactiva generada en la subestación
var Pgd{Ob}>= 0; 	   # Variable que indica que potencia activa va a ser inyectada por la GD
var Qgd{Ob};		   # Variable que indica que potencia reactiva va a ser inyectada por la GD

# var Wgd{Ob} binary;
#---------------------------------------------------------------------------
# Función objetivo
#----------------------------------------------------------------------------
# Minimizar pérdidas
#minimize FuncionObjetivo:
#    sum{ (i,j) in Ol }( R[i,j] * I[i,j]^2 );
# Maximizar Generación distribuida 
maximize FuncionObjetivo:
	sum{ i in Ob }(Pgd[i] );
#----------------------------------------------------------------------------
# Restricciones
#----------------------------------------------------------------------------
# Balance de potencia activa
subject to BalancePotenciaActiva{ i in Ob }:
	sum{ (k,i) in Ol }(P[k,i]) - sum{ (i,j) in Ol }( P[i,j] + R[i,j]*I[i,j]^2 ) + Ps[i] + Pgd[i] = PD[i];
# Balance de potencia reactiva
subject to BalancePotenciaReactiva{ i in Ob }:
	sum{ (k,i) in Ol }(Q[k,i]) - sum{ (i,j) in Ol }( Q[i,j] + X[i,j]*I[i,j]^2 ) + Qs[i] + Qgd[i] = QD[i];
# Caida de tensión en las líneas
subject to CaidaTension { (i,j) in Ol }:
	V[i]^2 - 2 * ( R[i,j] * P[i,j] + X[i,j] * Q[i,j] ) - (R[i,j]^2	+ X[i,j]^2) * I[i,j]^2 - V[j]^2  = 0;
# Potencia aparante
subject to PotenciaAparente{ (i,j) in Ol}:
	I[i,j]^2 * V[j]^2 = P[i,j]^2 + Q[i,j]^2;
# Límite de corriente en las líneas
subject to LimiteCorriente{ (i,j) in Ol }:
	0 <= I[i,j] <= Imax[i,j];
# Límite de la tensión
subject to LimiteTension{ i in Ob }:
	Vmin <= V[i] <= Vmax;
#----------------------------------------------------------------------------
# Restricciones relacionadas al planeamiento de GDs
#----------------------------------------------------------------------------
subject to generaciondistribuida1 { i in Ob}:
	Pgd[i] <= (Ongd[i]) * Pmax_gd; # Restricción para garantizar que los generadores no sean mayores a 1 MW
# subject to generaciondistribuida2:
subject to generaciondistribuida2 { i in Ob}:
	Qgd[i] = Pgd[i] * tan(acos(fp));
# Restricciones relacionadas al límite de pérdidas reconocidas
# subject to LimitePerdidas:
#   sum{(i,j) in Ol} (R[i,j] * I[i,j]^2) <= Ploss_max; # Restricción para garantizar que las pérdidas no superen el límite
# Restricción del desvío de tensión por desconexión abrupta del GD
#subject to DesvioTensionGD {i in Ob}:
#   V[i] >= Vmin_GD + Wgd[i]  * (Vnom - Vmin_GD);