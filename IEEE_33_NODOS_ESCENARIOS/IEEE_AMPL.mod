# IEEE_AMPL.mod - Modelo de optimización para generación solar con variabilidad horaria

# Definir el conjunto de horas del día
set H := 5..19;  # Definir el conjunto de horas correctamente


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
param Imax{Ol};		         # Corriente máxima de la línea
param linea{Ol};	         # Número de la línea
param PD{Ob};		         # Potencia activa demandada
param QD{Ob};		         # Potencia reactiva demandada
param QC{Ob};		         # Potencia reactiva del condensador en la barra
param Ongd{Ob};		         # Indica si hay generación distribuida en la barra
#---------------------------------------------------------------------------
# Parámetros variables por hora (dependientes del tiempo)
#---------------------------------------------------------------------------
param Demanda {Ob, H};       # Demanda por nodo y hora
param Generacion {Ob, H};    # Generación solar por nodo y hora
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
#---------------------------------------------------------------------------
# Variables del sistema
#---------------------------------------------------------------------------
var I{Ol, H};			   # Corriente en las líneas
var P{Ol, H};		       # Potencia activa en las líneas
var Q{Ol, H};		       # Potencia reactiva en las líneas
var V{Ob, H};			   # Tensión en las barras
var Ps{Ob, H} >= 0;	       # Potencia activa generada en la subestación
var Qs{Ob, H} >= 0;	       # Potencia reactiva generada en la subestación
var Pgd{Ob}>= 0; 	       # Variable que indica que potencia activa va a ser inyectada por la GD
var Qgd{Ob};		       # Variable que indica que potencia reactiva va a ser inyectada por la GD
#---------------------------------------------------------------------------
# Función objetivo
maximize FuncionObjetivo: 
	sum{ i in Ob, h in H}(Pgd[i]*Generacion[i,h] ); # Maximizar la energía generada por GDs considerando el perfil de generación
#----------------------------------------------------------------------------
# Restricciones
#----------------------------------------------------------------------------
# Balance de potencia activa
subject to BalancePotenciaActiva{ i in Ob, h in H }:
	sum{ (k,i) in Ol}( P[k,i,h] ) - sum{ (i,j) in Ol}( P[i,j,h] + R[i,j] * I[i,j,h]^2 ) + Ps[i,h] + Pgd[i]*Generacion[i,h] = PD[i]*Demanda[i,h];
# Balance de potencia reactiva
subject to BalancePotenciaReactiva{ i in Ob , h in H}:
	sum{ (k,i) in Ol}( Q[k,i,h] ) - sum{ (i,j) in Ol}( Q[i,j,h] + X[i,j] * I[i,j,h]^2 ) + Qs[i,h] + Qgd[i]*Generacion[i,h] = QD[i]*Demanda[i,h];
# Caida de tensión en las líneas
subject to CaidaTension { (i,j) in Ol, h in H}:
	V[i,h]^2 - 2 * ( R[i,j] * P[i,j,h] + X[i,j] * Q[i,j,h] ) - (R[i,j]^2	+ X[i,j]^2) * I[i,j,h]^2 - V[j,h]^2  = 0;
# Potencia aparante
subject to PotenciaAparente{ (i,j) in Ol, h in H}:
	I[i,j,h]^2 * V[j,h]^2 = P[i,j,h]^2 + Q[i,j,h]^2;
# Límite de corriente en las líneas
subject to LimiteCorriente{ (i,j) in Ol, h in H }:
	0 <= I[i,j,h] <= Imax[i,j];
# Límite de la tensión
subject to LimiteTension{ i in Ob, h in H }:
	Vmin <= V[i,h] <= Vmax;
# Fijar la potencia activa y reactiva en todos los nodos excepto en la subestación
subject to Fijar_Ps_Qs {i in Ob, h in H: i != Barra_SE}:
    Ps[i,h] = 0;
subject to Fijar_Qs {i in Ob, h in H: i != Barra_SE}:
    Qs[i,h] = 0;
# Fijar el voltaje de la subestación en Vnom para todas las horas
subject to Fijar_Voltaje_SE {h in H}:
    V[Barra_SE,h] = Vnom;
#---------------------------------------------------------------------------
# Restricciones relacionadas al planeamiento de generación distribuida (GD)
#---------------------------------------------------------------------------
subject to generaciondistribuida1 { i in Ob}:
	Pgd[i] <= (Ongd[i]) * Pmax_gd;# Restricción para garantizar que los GDs generen dentro de los límites establecidos
subject to generaciondistribuida2 { i in Ob}: # Relación entre potencia reactiva y activa según el factor de potencia
	Qgd[i] = Pgd[i]*tan(acos(fp));
