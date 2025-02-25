#---------------------------------------------------------------------------
# Flujo de potencia no lineal formulado como un problema de
# optimización matemática
#---------------------------------------------------------------------------
# definir los conjuntos con el comando: set
#---------------------------------------------------------------------------
set Ob;						 # conjunto de barras
set Ol within Ob cross Ob;	 # conjunto de lineas del sistema
#---------------------------------------------------------------------------
# Parametros
#---------------------------------------------------------------------------
param R{Ol};		         # Resistencia de la línea
param X{Ol};		         # Reactancia de la línea
param Z2{Ol};		         # Magnitud cuadrada de la impedancia de la línea
param G{Ol};		         # Conductancia de la línea
param B{Ol};		         # Subceptancia de la línea
param Imax{Ol};		         # Corriente máxima de la línea
param linea{Ol};	         # Número de la línea
param PD{Ob};		         # Potencia activa demandada
param QD{Ob};		         # Potencia reactiva demandada
param QC{Ob};		         # Potencia reactiva del condensador en la barra
param Ongd{Ob};		         # Indica si hay generación distribuida en la barra

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
var L{Ol}  >= 0;	         # Magnitud de Corriente cuadrada en las líneas
var P{Ol};		             # Potencia activa en las líneas
var Q{Ol};		             # Potencia reactiva en las líneas
var U{Ob}  >= 0;             # Magnitud de Tensión cuadrada en las barras 
var Ps{Ob} >= 0;	         # Potencia activa generada en la subestación
var Qs{Ob} >= 0;	         # Potencia reactiva generada en la subestación
var Pgd{Ob}>= 0; 	         # Variable que indica que potencia activa va a ser inyectada por la GD
var Qgd{Ob};		   # Variable que indica que potencia reactiva va a ser inyectada por la GD
#---------------------------------------------------------------------------
# Función objetivo
#----------------------------------------------------------------------------
# Maximizar Generación distribuida 
maximize FuncionObjetivo:
	sum{ i in Ob }(Pgd[i] );
#----------------------------------------------------------------------------
# Restricciones
#----------------------------------------------------------------------------
# Balance de potencia activa
subject to BalancePotenciaActiva{ i in Ob }:
    sum{ (k,i) in Ol }(P[k,i]) - sum{ (i,j) in Ol }( P[i,j] + R[i,j]*L[i,j] ) + Ps[i] + Pgd[i] = PD[i];
# Balance de potencia reactiva
subject to BalancePotenciaReactiva{ i in Ob }:
	sum{ (k,i) in Ol }(Q[k,i]) - sum{ (i,j) in Ol }( Q[i,j] + X[i,j]*L[i,j] ) + Qs[i]  + Qgd[i]= QD[i];
# Caida de tensión en las líneas
subject to CaidaTension { (i,j) in Ol }:
    U[i] - 2 * ( R[i,j] * P[i,j] + X[i,j] * Q[i,j] ) - (R[i,j]^2 + X[i,j]^2) * L[i,j] - U[j] = 0;
# Límite de corriente en las líneas
subject to LimiteCorriente{ (i,j) in Ol }:
    L[i,j] <= Imax[i,j]^2;
# Límite de la tensión
subject to LimiteTension{ i in Ob }:
	Vmin^2 <= U[i] <= Vmax^2;
# Fijar la potencia activa y reactiva en todos los nodos excepto en la subestación
subject to Fijar_Ps_Qs {i in Ob: i != Barra_SE}:
    Ps[i] = 0;
subject to Fijar_Qs {i in Ob: i != Barra_SE}:
    Qs[i] = 0;
# Fijar el voltaje de la subestación en Vnom para todas las horas
subject to Fijar_Voltaje_SE :
    U[Barra_SE] = Vnom^2;
#----------------------------------------------------------------------------
# Restricciones relacionadas al planeamiento de GDs
#----------------------------------------------------------------------------
subject to generaciondistribuida1 { i in Ob}:
	Pgd[i] <= (Ongd[i]) * Pmax_gd; # Restricción para garantizar que los generadores no sean mayores a 1 MW
	# subject to generaciondistribuida2:
subject to generaciondistribuida2 { i in Ob}:
	Qgd[i] = Pgd[i] * tan(acos(fp));
#----------------------------------------------------------------------------
# Restricciones relacionadas con la relajación cónica convexa de segundo orden
#----------------------------------------------------------------------------
subject to relajacionconica { (i,j) in Ol }:
    (2 * P[i,j])^2 + (2 * Q[i,j])^2 + (U[i] - L[i,j])^2 <= (U[i] + L[i,j])^2;

	