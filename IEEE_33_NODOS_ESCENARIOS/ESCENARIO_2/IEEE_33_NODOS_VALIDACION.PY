import os
import pandas as pd
from amplpy import AMPL, ampl_notebook

os.system('cls' if os.name == 'nt' else 'clear')
main_folder = "AMPL"
sub_folder = "IEEE_33_NODOS_ESCENARIOS"
sub_sub_folder = "ESCENARIO_2"
base_path = os.getcwd()
target_path = os.path.join(base_path, main_folder, sub_folder, sub_sub_folder)
os.chdir(target_path)

print("Directorio actual:", os.getcwd())
# Inicializar AMPL
ampl = ampl_notebook(license_uuid="0b2bd731-acd3-469f-9b23-57614ace54e3")


def determinar_peores_horas():
    demanda = pd.read_csv('IEEE_AMPL_PERFILES_DEMANDA.dat', sep=' ', skiprows=1, header=None, names=["Nodo", "Hora", "Valor"])
    generacion = pd.read_csv('IEEE_AMPL_PERFILES_GENERACION.dat', sep=' ', skiprows=1, header=None, names=["Nodo", "Hora", "Valor"])
    demanda = demanda.dropna().reset_index(drop=True)
    generacion = generacion.dropna().reset_index(drop=True)
    demanda = demanda[demanda["Nodo"] != ";"]
    generacion = generacion[generacion["Nodo"] != ";"]
    demanda["Nodo"]     = pd.to_numeric(demanda["Nodo"], errors='coerce')
    demanda["Hora"]     = pd.to_numeric(demanda["Hora"], errors='coerce')
    demanda["Valor"]    = pd.to_numeric(demanda["Valor"], errors='coerce')
    generacion["Nodo"]  = pd.to_numeric(generacion["Nodo"], errors='coerce')
    generacion["Hora"]  = pd.to_numeric(generacion["Hora"], errors='coerce')
    generacion["Valor"] = pd.to_numeric(generacion["Valor"], errors='coerce')
    demanda = demanda.dropna().reset_index(drop=True)
    generacion = generacion.dropna().reset_index(drop=True)
    demanda = demanda.drop_duplicates(subset=["Nodo", "Hora"])
    generacion = generacion.drop_duplicates(subset=["Nodo", "Hora"])
    demanda_pivot = demanda.pivot(index="Nodo", columns="Hora", values="Valor").fillna(0)
    generacion_pivot = generacion.pivot(index="Nodo", columns="Hora", values="Valor").fillna(0)
    #print("✅ Demanda Pivotada:\n", demanda_pivot)
    #print("✅ Generación Pivotada:\n", generacion_pivot)
    balance_hora = {}
    for h in range(1, 24):
        if h in generacion_pivot.columns and h in demanda_pivot.columns:
            balance_hora[h] = generacion_pivot[h].sum() - (demanda_pivot[h] * demanda_pivot[1]).sum()
        else:
            print(f"⚠️ La hora {h} no está en los datos de generación o demanda")
    peores_horas = sorted(balance_hora, key=balance_hora.get, reverse=True)[:6]
    print(f"🔥 Peores horas seleccionadas: {peores_horas}")
    return peores_horas

peores_horas = determinar_peores_horas()
horas_str = " ".join(map(str, peores_horas))
with open('IEEE_AMPL_PEORES_HORAS.dat', 'w') as f:
    f.write(f"set PH := {horas_str};\n")


def solve_optimization():
    ampl = AMPL()
    ampl.reset()

    # Cargar el modelo y los datos
    ampl.read('IEEE_AMPL.mod')  # Primero el modelo
    ampl.readData('IEEE_AMPL.dat')  # Datos del sistema (líneas, barras)
    ampl.readData('IEEE_AMPL_PERFILES_GENERACION.dat')  # Perfiles de demanda y generación
    ampl.readData('IEEE_AMPL_PERFILES_DEMANDA.dat')  # Perfiles de demanda y generación
    ampl.readData('IEEE_AMPL_PEORES_HORAS.dat')  # Peores horas

    ampl.display("Ongd")
    ampl.display("Pmax_gd")
    ampl.display("H")

    # Configurar el solver (BARON, CONOPT, CPLEX, GUROBI, XPRESS, IPOPT, COPT & MOSEK)  Couenne
    ampl.setOption("solver", "ipopt")
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
          f.write(f"New Generator.{bus_number} Bus1={bus_number} Phases=3 kV=12.66 kW={pgd_value:.3f} kvar=0\n")
        print("✅ Archivo DSS generado con éxito.")
        # Calcular HC
        HC = pgd_values["Pgd.val"].sum()  # Suma de los valores de Pgd
        print(f"✅ El valor de HC es: {HC}")
    except RuntimeError:
        print("⚠️ No se pudieron obtener valores para la variable 'Pgd'.")

    try:
       # Obtener los valores de V desde AMPL
       v_values = ampl.getVariable("V").getValues().toPandas()
       v_values.reset_index(inplace=True)
       print("📊 Valores de la variable v:")
       print(v_values)
       # Obtener el mínimo y el máximo en pu
       v_min_pu = v_values["V.val"].min()
       v_max_pu = v_values["V.val"].max()
       # Imprimir los valores correctos
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
            values.reset_index(inplace=True)
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
        "V_max": v_max_pu
    }
if __name__ == "__main__":
    resultados = solve_optimization()
    # print("Resultados del problema de optimización:", resultados)

