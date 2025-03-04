import numpy as np
import pandas as pd
import random
import os
import matplotlib.pyplot as plt

os.system('cls' if os.name == 'nt' else 'clear')
main_folder = "AMPL"
sub_folder = "IEEE_33_NODOS_ESCENARIOS"
sub_sub_folder = "ESCENARIO_4"
base_path = os.getcwd()
target_path = os.path.join(base_path, main_folder, sub_folder, sub_sub_folder)
os.chdir(target_path)
print("📂 Directorio actual:", os.getcwd())


# Crear directorio para gráficos si no existe
graphs_path = os.path.join(target_path, "GRÁFICOS")
os.makedirs(graphs_path, exist_ok=True)

perfiles_demanda = {
    "Demanda_Promedio_1": [0.29, 0.31, 0.26, 0.27, 0.31, 0.33, 0.38, 0.37, 0.41, 0.47, 0.53, 0.55, 0.54, 0.50, 0.47, 0.42, 0.52, 0.61, 0.92, 1.00, 0.90, 0.75, 0.63, 0.45,1],
    "Demanda_Promedio_2": [0.30, 0.27, 0.26, 0.27, 0.24, 0.41, 0.44, 0.36, 0.36, 0.45, 0.41, 0.45, 0.56, 0.39, 0.41, 0.45, 0.47, 0.52, 0.77, 0.91, 0.98, 1.00, 0.64, 0.37,1],
    "Demanda_Promedio_3": [0.28, 0.25, 0.26, 0.25, 0.23, 0.44, 0.69, 0.44, 0.44, 0.36, 0.38, 0.51, 0.40, 0.38, 0.37, 0.41, 0.49, 0.34, 0.61, 0.81, 1.00, 0.86, 0.69, 0.36,1],
    "Demanda_Promedio_4": [0.45, 0.41, 0.38, 0.44, 0.39, 0.52, 0.68, 0.68, 0.59, 0.59, 0.63, 0.54, 0.59, 0.62, 0.60, 0.61, 0.48, 0.54, 0.72, 0.82, 1.00, 0.78, 0.67, 0.58,1],
    "Demanda_Promedio_5": [0.33, 0.38, 0.38, 0.34, 0.34, 0.36, 0.37, 0.57, 0.62, 0.52, 0.56, 1.00, 0.95, 0.74, 0.66, 0.65, 0.54, 0.48, 0.62, 0.73, 0.83, 0.73, 0.67, 0.55,1],
    #"Demanda_Promedio_6": [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.34, 0.76, 0.91, 0.93, 0.72, 0.41, 0.70, 0.95, 1.00, 0.98, 0.85, 0.61, 0.16, 0.07, 0.06, 0.05,1]
}
perfiles_generacion = {
    "Generacion_Tipica_1": [0, 0, 0, 0, 0, 0.10, 0.30, 0.50, 0.70, 0.85, 0.95, 1.00, 0.98, 0.92, 0.85, 0.70, 0.50, 0.30, 0.10,  0, 0, 0, 0, 0,1],  # Día despejado
    #"Generacion_Tipica_2": [0, 0, 0, 0, 0, 0.05, 0.20, 0.40, 0.60, 0.75, 0.85, 0.88, 0.85, 0.75, 0.60, 0.40, 0.20, 0.05,    0,  0, 0, 0, 0, 0,1],  # Día nublado
    #"Generacion_Tipica_3": [0, 0, 0, 0, 0, 0.08, 0.25, 0.45, 0.65, 0.82, 0.90, 0.94, 0.92, 0.88, 0.78, 0.62, 0.45, 0.20, 0.08,  0, 0, 0, 0, 0,1]   # Día con nubosidad parcial
}

# Definir parámetros
num_dias = 30
num_horas_dia = 24
num_nodos = 33  # Número de nodos en el sistema

# Almacenar datos de demanda
demanda_por_nodo = {nodo: [] for nodo in range(2, num_nodos + 1)}

# Escribir archivos DAT con variaciones diarias
with open("IEEE_AMPL_PERFILES_GENERACION.dat", "w") as f_gen, open("IEEE_AMPL_PERFILES_DEMANDA.dat", "w") as f_dem:
    f_gen.write("\nparam Generacion :=\n")
    f_dem.write("\nparam Demanda :=\n")
    
    for dia in range(1, num_dias + 1):
        for nodo in range(2, num_nodos + 1):
            perfil_demanda = random.choice(list(perfiles_demanda.values()))
            perfil_generacion = random.choice(list(perfiles_generacion.values()))
            
            for hora in range(1, num_horas_dia + 1):
                variacion_demanda = random.uniform(0.9, 1.1)
                variacion_generacion = random.uniform(0.9, 1.1)
                
                demanda_valor = round(perfil_demanda[hora - 1] * variacion_demanda, 4)
                generacion_valor = round(perfil_generacion[hora - 1] * variacion_generacion, 4)
                
                f_dem.write(f"{dia} {nodo} {hora} {demanda_valor}\n")
                f_gen.write(f"{dia} {nodo} {hora} {generacion_valor}\n")
                
                # Almacenar demanda para graficar
                demanda_por_nodo[nodo].append(demanda_valor)
    
    f_gen.write(";\n")
    f_dem.write(";\n")

# Generar y guardar gráficos de demanda para cada nodo
for nodo, demanda in demanda_por_nodo.items():
    plt.figure(figsize=(12, 6))
    plt.plot(range(1, num_dias * num_horas_dia + 1), demanda, label=f"Nodo {nodo}", color='b')
    plt.xlabel("Horas en el Año")
    plt.ylabel("Demanda (p.u.)")
    plt.title(f"Perfil de Demanda del Nodo {nodo} a lo Largo del Año")
    plt.legend()
    plt.grid(True)
    
    # Guardar la imagen
    plt.savefig(os.path.join(graphs_path, f"Demanda_Nodo_{nodo}.png"), dpi=300)
    plt.close()

print("\n✅ Archivos IEEE_AMPL_PERFILES_GENERACION.dat e IEEE_AMPL_PERFILES_DEMANDA.dat generados correctamente para 360 días con variaciones diarias.")
print("📊 Gráficos de demanda guardados en la carpeta 'graficos'.")