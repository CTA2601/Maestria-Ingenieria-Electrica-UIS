import random
import pandas as pd

perfiles_demanda = {
    #"Demanda_Promedio_1": [0.29, 0.31, 0.26, 0.27, 0.31, 0.33, 0.38, 0.37, 0.41, 0.47, 0.53, 0.55, 0.54, 0.50, 0.47, 0.42, 0.52, 0.61, 0.92, 1.00, 0.90, 0.75, 0.63, 0.45,1],
    #"Demanda_Promedio_2": [0.30, 0.27, 0.26, 0.27, 0.24, 0.41, 0.44, 0.36, 0.36, 0.45, 0.41, 0.45, 0.56, 0.39, 0.41, 0.45, 0.47, 0.52, 0.77, 0.91, 0.98, 1.00, 0.64, 0.37,1],
    #"Demanda_Promedio_3": [0.28, 0.25, 0.26, 0.25, 0.23, 0.44, 0.69, 0.44, 0.44, 0.36, 0.38, 0.51, 0.40, 0.38, 0.37, 0.41, 0.49, 0.34, 0.61, 0.81, 1.00, 0.86, 0.69, 0.36,1],
    #"Demanda_Promedio_4": [0.45, 0.41, 0.38, 0.44, 0.39, 0.52, 0.68, 0.68, 0.59, 0.59, 0.63, 0.54, 0.59, 0.62, 0.60, 0.61, 0.48, 0.54, 0.72, 0.82, 1.00, 0.78, 0.67, 0.58,1],
    "Demanda_Promedio_5": [0.33, 0.38, 0.38, 0.34, 0.34, 0.36, 0.37, 0.57, 0.62, 0.52, 0.56, 1.00, 0.95, 0.74, 0.66, 0.65, 0.54, 0.48, 0.62, 0.73, 0.83, 0.73, 0.67, 0.55,1],
    #"Demanda_Promedio_6": [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.34, 0.76, 0.91, 0.93, 0.72, 0.41, 0.70, 0.95, 1.00, 0.98, 0.85, 0.61, 0.16, 0.07, 0.06, 0.05,1]
}
perfiles_generacion = {
    "Generacion_Tipica_1": [0, 0, 0, 0, 0, 0.10, 0.30, 0.50, 0.70, 0.85, 0.95, 1.00, 0.98, 0.92, 0.85, 0.70, 0.50, 0.30, 0.10,  0, 0, 0, 0, 0,1],  # Día despejado
    #"Generacion_Tipica_2": [0, 0, 0, 0, 0, 0.05, 0.20, 0.40, 0.60, 0.75, 0.85, 0.88, 0.85, 0.75, 0.60, 0.40, 0.20, 0.05,    0,  0, 0, 0, 0, 0,1],  # Día nublado
    #"Generacion_Tipica_3": [0, 0, 0, 0, 0, 0.08, 0.25, 0.45, 0.65, 0.82, 0.90, 0.94, 0.92, 0.88, 0.78, 0.62, 0.45, 0.20, 0.08,  0, 0, 0, 0, 0,1]   # Día con nubosidad parcial
}

# Número de nodos
N = 33
H = 25

# Asignar perfiles aleatorios a cada nodo
asignacion_demanda = {nodo: random.choice(list(perfiles_demanda.items())) for nodo in range(1, N+1)}
asignacion_generacion = {nodo: random.choice(list(perfiles_generacion.items())) for nodo in range(1, N+1)}

# Escribir en IEEE_AMPL.dat
with open("IEEE_AMPL_PERFILES.dat", "w") as f:
    # Escribir perfiles de generación
    f.write("\nparam Generacion :=\n")
    for nodo, (perfil_nombre, perfil_valores) in asignacion_generacion.items():
        for hora in range(1, H+1):
            f.write(f"{nodo} {hora} {perfil_valores[hora-1]}\n")
    f.write(";\n")
    
    # Escribir perfiles de demanda
    f.write("\nparam Demanda :=\n")
    for nodo, (perfil_nombre, perfil_valores) in asignacion_demanda.items():
        for hora in range(1, H+1):
            f.write(f"{nodo} {hora} {perfil_valores[hora-1]}\n")
    f.write(";\n")


# Crear DataFrames para mostrar las asignaciones
df_demanda = pd.DataFrame([(nodo, perfil_nombre) for nodo, (perfil_nombre, _) in asignacion_demanda.items()], columns=["Nodo", "Perfil Demanda"])
df_generacion = pd.DataFrame([(nodo, perfil_nombre) for nodo, (perfil_nombre, _) in asignacion_generacion.items()], columns=["Nodo", "Perfil Generación"])

# Imprimir tablas de asignación
print("\n📊 Asignación de Perfiles de Demanda:")
print(df_demanda)
print("\n📊 Asignación de Perfiles de Generación:")
print(df_generacion)

print("\n✅ Archivo IEEE_AMPL.dat generado con asignación aleatoria de perfiles de demanda y generación.")