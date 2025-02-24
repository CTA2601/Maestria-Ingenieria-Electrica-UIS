import os
import pandas as pd
import numpy as np

os.system('cls' if os.name == 'nt' else 'clear')
main_folder = "AMPL"
sub_folder = "IEEE_33_NODOS"
sub_sub_folder = "SOCP"
base_path = os.getcwd()
target_path = os.path.join(base_path, main_folder, sub_folder, sub_sub_folder)
os.chdir(target_path)
print("📂 Directorio actual:", os.getcwd())

from amplpy import AMPL, ampl_notebook
ampl = ampl_notebook(
    license_uuid="0b2bd731-acd3-469f-9b23-57614ace54e3") # your license UUID (e.g., free ampl.com/ce or ampl.com/courses licenses)


def solve_optimization():
    ampl = AMPL()
    ampl.reset()

    ampl.read('IEEE_AMPL.mod') 
    ampl.readData('IEEE_AMPL.dat')

    # Configurar el solver (BARON, CONOPT, CPLEX, GUROBI, XPRESS, IPOPT, COPT & MOSEK)  Couenne
    ampl.setOption("solver", "cplex")
    ampl.solve("FuncionObjetivo")
    try:
       # Obtener los valores de V desde AMPL
       U_values = ampl.getVariable("U").getValues().toPandas()
       U_values.reset_index(inplace=True)
       print("📊 Valores de la variable U:")
       print(U_values)
       # Obtener el mínimo y el máximo en pu
       V_values = np.sqrt(U_values["U.val"] * np.sqrt(2))
       print(V_values)
       v_min_pu = V_values.min()
       v_max_pu = V_values.max()
       # Imprimir los valores correctos
       print(f"⚡ Tensión mínima (pu): {v_min_pu:.4f}")
       print(f"⚡ Tensión máxima (pu): {v_max_pu:.4f}")
    except RuntimeError:
       print("⚠️ No se pudieron obtener valores para la variable 'V'.")
  
    # Obtener el valor de la función objetivo
    funcion_objetivo = ampl.getObjective("FuncionObjetivo").value()

    # Retornar resultados
    return {
        "FuncionObjetivo": funcion_objetivo,
        "V_min": v_min_pu,
        "V_max": v_max_pu
    }

if __name__ == "__main__":
    resultados = solve_optimization()
    # print("Resultados del problema de optimización:", resultados)
