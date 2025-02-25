import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from amplpy import AMPL, ampl_notebook

# Limpiar pantalla
os.system('cls' if os.name == 'nt' else 'clear')

# Definir rutas
directorio_principal = "AMPL"
directorio_sub = "IEEE_33_NODOS"
directorio_subsub = "SOCP"
base_path = os.getcwd()
directorio_destino = os.path.join(base_path, directorio_principal, directorio_sub, directorio_subsub)
os.chdir(directorio_destino)
print("📂 Directorio actual:", os.getcwd())

# Inicializar AMPL
ampl = ampl_notebook(license_uuid="0b2bd731-acd3-469f-9b23-57614ace54e3")

def solve_optimization():
    ampl = AMPL()
    ampl.reset()

    # Cargar el modelo y los datos
    ampl.read('IEEE_AMPL.mod')  # Primero el modelo
    ampl.readData('IEEE_AMPL.dat')  # Datos del sistema (líneas, barras)


    #ampl.display("Ongd")
    ampl.display("Pmax_gd")

    # Configurar el solver (BARON, CONOPT, CPLEX, GUROBI, XPRESS, IPOPT, COPT & MOSEK)  Couenne
    ampl.option['solver'] = 'ipopt'  # Mejor para no lineales
    ampl.solve("FuncionObjetivo")
    # Obtener valor de Pgd
    try:
        pgd_values = ampl.getVariable("Pgd").getValues().toPandas()
        pgd_values.reset_index(inplace=True)
        print("📊 Valores de la variable Pgd:")
        print(pgd_values)

        # Crear el archivo DSS con los valores de Pgd
        with open("IEEE_33_GENERADORES.dss", "w") as f:
         for _, row in pgd_values.iterrows():
          bus_number = int(row['index']) 
          pgd_value = row['Pgd.val'] * 1000
          f.write(f"New Generator.{bus_number} Bus1={bus_number} Phases=3 kV=12.66 kW={pgd_value:.3f} kvar=0  Vminpu=0.7 Vmaxpu=1.2 conn=Wye Model=1\n")
        print("✅ Archivo DSS generado con éxito.")
        # Calcular HC
        HC = pgd_values["Pgd.val"].sum()  # Suma de los valores de Pgd
        print(f"✅ El valor de HC es: {HC}")
    except RuntimeError:
        print("⚠️ No se pudieron obtener valores para la variable 'Pgd'.")

    try:
       # Obtener los valores de V desde AMPL
        U_values = ampl.getVariable("U").getValues().toPandas()
        U_values.reset_index(inplace=True)
        print("📊 Valores de la variable U:")
        #print(U_values)
        # Cálculo de tensiones en p.u.
        V_values =np.sqrt(U_values["U.val"])
        print(V_values)
        v_min_pu = V_values.min()
        v_max_pu = V_values.max()
        print(f"⚡ Tensión mínima (pu): {v_min_pu:.4f}")
        print(f"⚡ Tensión máxima (pu): {v_max_pu:.4f}")
    except RuntimeError:
       print("⚠️ No se pudieron obtener valores para la variable 'V'.")
  
    # Obtener el valor de la función objetivo
    funcion_objetivo = ampl.getObjective("FuncionObjetivo").value()

    # Obtener y mostrar los valores de las variables de decisión
    variables = ampl.getVariables()  # Devuelve un EntityMap
    all_values = []
    for var_name, variable in variables:  # Desempaqueta las tuplas (nombre, objeto_variable)
        try:
            # Obtener los valores de la variable y convertirlos a un DataFrame
            values = variable.getValues().toPandas()
            values.reset_index(inplace=True)        # Mantener todos los índices
            values.insert(0, "Variable", var_name)  # Agregar columna Variable
            all_values.append(values)
        except RuntimeError:
            print(f"⚠️ No se pudieron obtener valores para la variable '{var_name}'.")

    # Combinar todos los valores en un único DataFrame y exportar a CSV
    if all_values:
        combined_values = pd.concat(all_values, ignore_index=True)
        combined_values.to_csv("IEEE_33_SOLVGD.csv", index=False)
        print("✅ Los valores de las variables se han exportado a 'IEEE_33_SOLVGD.csv'.")
    else:
        print("📢 No se encontraron valores de variables para exportar.")

    # Retornar resultados
    return {
        "FuncionObjetivo": funcion_objetivo,
        "HC": HC,
        "V_min": v_min_pu,
        "V_max": v_max_pu,
        "V": V_values.tolist()
    }

if __name__ == "__main__":
    resultados = solve_optimization()
    # print("Resultados del problema de optimización:", resultados)

# Graficar voltajes en los nodos
if resultados:
    
    # Obtener nodos con voltaje máximo y mínimo
    max_voltage = max(resultados["V"])
    min_voltage = min(resultados["V"])
    max_bus = resultados["V"].index(max_voltage) + 1
    min_bus = resultados["V"].index(min_voltage) + 1
    
    # Imprimir resultados
    print(f"⚡ Voltaje Máximo: {max_voltage:.4f} p.u. en el nodo {max_bus}")
    print(f"⚡ Voltaje Mínimo: {min_voltage:.4f} p.u. en el nodo {min_bus}")
    
    
    plt.figure(figsize=(12, 6))
    plt.plot(range(1, len(resultados["V"]) + 1), resultados["V"], marker='o', linestyle='-', color='b', label="Voltaje en p.u.")
    plt.axhline(y=1.05, color='r', linestyle='dashed', label="Límite Superior (1.05 p.u.)")
    plt.axhline(y=0.95, color='b', linestyle='dashed', label="Límite Inferior (0.95 p.u.)")
    plt.ylim(0.90, 1.10)
    plt.xlabel("Nodos del Sistema", fontsize=12)
    plt.ylabel("Voltaje (p.u.)", fontsize=12)
    plt.title("Tensión en todos los nodos del sistema", fontsize=14)
    plt.xticks(range(1, len(resultados["V"]) + 1), rotation=90, fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.show()

