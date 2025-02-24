# IEEE_AMPL.mod - Modelo de optimización para generación solar con variabilidad horaria

# Definir el conjunto de horas del día
set PH;  # Definir el conjunto de horas correctamente


# Parámetros
#---------------------------------------------------------------------------
# Flujo de potencia no lineal formulado como un problema de
# optimización matemática
#---------------------------------------------------------------------------
# definir los conjuntos con el comando: set
#---------------------------------------------------------------------------
set Ob;						           # conjunto de barras
set Ol within Ob cross Ob;	           # conjunto de lineas del sistema
#---------------------------------------------------------------------------
# Parámetros variables por nodo
#---------------------------------------------------------------------------
param R{Ol};		         # Resistencia de la línea
param X{Ol};		         # Reactancia de la línea
param Z2{Ol};		         # Z^2 = R^2 + X^2
param G{Ol};		         # Resistencia de la línea
param B{Ol};		         # Reactancia de la línea
param Imax{Ol};		         # Corriente máxima de la línea
param linea{Ol};	         # Número de la línea
param PD{Ob};		         # Potencia activa demandada
param QD{Ob};		         # Potencia reactiva demandada
param QC{Ob};		         # Potencia reactiva del condensador en la barra
param Ongd{Ob};		         # Indica si hay generación distribuida en la barra
#---------------------------------------------------------------------------
# Parámetros variables por hora (dependientes del tiempo)
#---------------------------------------------------------------------------
param Demanda {Ob,PH};       # Demanda por nodo y hora
param Generacion {Ob,PH};    # Generación solar por nodo y hora
#---------------------------------------------------------------------------
# Parámetros fijos
#---------------------------------------------------------------------------
param Vnom; 		         # Tensión nominal de la red
param Vmin;			         # Tensión mínima
param Vmax;			         # Tensión máxima
param Barra_SE;		         # Barra de la subestación
param Pmax_gd := 10 ;        # Potencia maxima que puede tener los GDs
param fp := 1;	             # Factor de potencia de la GD
param Nmax_gd := 10;         # Número máximo de GDs
param Ploss_max:= 0.191084;  # Límite máximo de pérdidas permitidas en la red
param factor_qgd := tan(acos(fp));  # Se precalcula fuera de las restricciones
#---------------------------------------------------------------------------
# Variables del sistema
#---------------------------------------------------------------------------
var Ps{Ob,PH} >= 0;	   # Potencia activa generada en la subestación
var Qs{Ob,PH} >= 0;	   # Potencia reactiva generada en la subestación
var U{Ob, PH};         # Cambio de variable
var PR{Ol, PH};        # Parte Real SOCP
var PI{Ol, PH};        # Parte imaginaria SOCP
var Pgd{Ob} >= 0; 	   # Variable que indica que potencia activa va a ser inyectada por la GD
var Qgd{Ob};		   # Variable que indica que potencia reactiva va a ser inyectada por la GD

#---------------------------------------------------------------------------
# Función objetivo
maximize FuncionObjetivo: 
	sum{ i in Ob}((Pgd[i])*Ongd[i]); # Maximizar la capacidad norminal  de los GDs
#----------------------------------------------------------------------------
# Restricciones
#----------------------------------------------------------------------------
# Balance de potencia activa
subject to BalancePotenciaActiva {i in Ob, h in PH}:
   sum{(i,j) in Ol} ((G[i,j] * PR[i,j,h]) - (B[i,j] * PI[i,j,h]) -(G[i,j] * U[i,h])) + Ps[i,h] + Pgd[i] * Generacion[i,h] = PD[i] * Demanda[i,h];
# Balance de potencia reactiva
subject to BalancePotenciaReactiva {i in Ob, h in PH}:
   sum{(i,j) in Ol} ((B[i,j] * PR[i,j,h]) - (G[i,j] * PI[i,j,h])- (B[i,j] * U[i,h])) + Qs[i,h]  = QD[i] * Demanda[i,h];
# Cambio de Variable
subject to CambioVariable{ (i,j) in Ol, h in PH}:
    PR[i,j,h]^2 + PI[i,j,h]^2 <= U[i,h] * U[j,h];
# Límite de la tensión
subject to LimiteTension{ i in Ob, h in PH }:
    Vmin^2 <= U[i,h] <= Vmax^2;
# Fijar la potencia activa y reactiva en todos los nodos excepto en la subestación
subject to Fijar_Ps {i in Ob, h in PH: i != Barra_SE}:
    Ps[i,h] = 0;
subject to Fijar_Qs {i in Ob, h in PH: i != Barra_SE}:
    Qs[i,h] = 0;
# Fijar el voltaje de la subestación en Vnom para todas las horas
subject to Fijar_Voltaje_SE {i in Ob: i == Barra_SE, h in PH}:
    U[Barra_SE,h] = Vnom^2;
#---------------------------------------------------------------------------
# Restricciones relacionadas al planeamiento de generación distribuida (GD)
#---------------------------------------------------------------------------
# Restricción para garantizar que los GDs generen dentro de los límites establecidos
subject to generaciondistribuida1 { i in Ob}:
	Pgd[i] <= (Ongd[i]) * Pmax_gd;
# Relación entre potencia reactiva y activa según el factor de potencia
