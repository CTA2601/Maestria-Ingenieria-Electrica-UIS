#---------------------------------------------------------------------------
# Flujo de potencia cónica formulado como un problema de
# optimización matemática
#---------------------------------------------------------------------------
# definir los conjuntos con el comando: set
#---------------------------------------------------------------------------
set Ob;						           # conjunto de barras
set Ol within Ob cross Ob;	           # conjunto de lineas del sistema
#---------------------------------------------------------------------------
param R{Ol};		         # Resistencia de la línea
param X{Ol};		         # Reactancia de la línea
param Z2{Ol};		         # Z^2 = R^2 + X^2
param G{Ol};                 # Conductancia definida en el .dat  # Conductancia
param B{Ol};                 # Susceptancia definida en el .dat # Susceptancia
param Imax{Ol};		         # Corriente máxima de la línea
param linea{Ol};	         # Número de la línea
param PD{Ob};		         # Potencia activa demandada
param QD{Ob};		         # Potencia reactiva demandada
#---------------------------------------------------------------------------
param Vnom; 		         # Tensión nominal de la red
param Barra_SE;		         # Barra de la subestación
#---------------------------------------------------------------------------
# Se definen las variables con el comando var
#---------------------------------------------------------------------------
var U{Ob}  ;        # Cambio de variable
var PR{Ol} ;        # Parte Real SOCP
var PI{Ol} ;        # Parte imaginaria SOCP
#---------------------------------------------------------------------------
# Función objetivo
#----------------------------------------------------------------------------
maximize FuncionObjetivo:
    sum{ (i,j) in Ol } (PR[i,j]);
#----------------------------------------------------------------------------

# Balance de potencia reactiva# Balance de potencia activa
subject to BalancePotenciaActiva {i in Ob diff {1}}:
    - sum{(i,j) in Ol} ( sqrt(2) * U[i] * G[i,j] - (G[i,j] * PR[i,j] - B[i,j] * PI[i,j]) ) 
    - sum{(j,i) in Ol} ( sqrt(2) * U[j] * G[j,i] - (G[j,i] * PR[j,i] - B[j,i] * PI[j,i]) ) 
    = PD[i];
subject to BalancePotenciaReactiva {i in Ob diff {1}}:
    - sum{(i,j) in Ol} ( sqrt(2) * U[i] * B[i,j] - (B[i,j] * PR[i,j] - G[i,j] * PI[i,j]) ) 
    - sum{(j,i) in Ol} ( sqrt(2) * U[j] * B[j,i] - (B[j,i] * PR[j,i] - G[j,i] * PI[j,i]) ) 
    = QD[i];
# Cambio de Variable (Restricción Cónica)
subject to CambioVariable { (i,j) in Ol }:
    PR[i,j]^2 + PI[i,j]^2 <= 2*U[i] * U[j];
# Fijar el voltaje de la subestación en Vnom
subject to Fijar_Voltaje_SE:
    U[Barra_SE] = Vnom^2/ sqrt(2);
# Las tensiones son validas solo como valores positivos
subject to Tension {i in Ob diff {1}}:
    U[i] >= 0;
# Restricción de simetría para la potencia activa
#subject to SimetriaPR { (i,j) in Ol }:
#    PR[i,j] = PR[j,i];
# Restricción de anti-simetría para la potencia reactiva
#subject to SimetriaPI { (i,j) in Ol }:
#    PI[i,j] = -PI[j,i];

    