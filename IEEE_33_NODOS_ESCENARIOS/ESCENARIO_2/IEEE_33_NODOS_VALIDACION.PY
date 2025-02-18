import os
import pandas as pd
from amplpy import AMPL, ampl_notebook
import dss
import matplotlib.pyplot as plt
import numpy as np

os.system('cls' if os.name == 'nt' else 'clear')
main_folder = "AMPL"
sub_folder = "IEEE_33_NODOS_ESCENARIOS"
sub_sub_folder = "ESCENARIO_2"
base_path = os.getcwd()
target_path = os.path.join(base_path, main_folder, sub_folder, sub_sub_folder)
os.chdir(target_path)

print("📂 Directorio actual:", os.getcwd())
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
    for h in range(1, 25):
        if h in generacion_pivot.columns and h in demanda_pivot.columns:
            balance_hora[h] = generacion_pivot[h].sum() - (demanda_pivot[h] * demanda_pivot[1]).sum()
        else:
            print(f"⚠️ La hora {h} no está en los datos de generación o demanda")
    peores_horas = sorted(balance_hora, key=balance_hora.get, reverse=True)[:1]
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

    #ampl.display("Ongd")
    #ampl.display("Pmax_gd")
    #ampl.display("PH")

    # Configurar el solver (BARON, CONOPT, CPLEX, GUROBI, XPRESS, IPOPT, COPT & MOSEK)  Couenne
    ampl.setOption("solver", "ipopt")
    ampl.solve("FuncionObjetivo")
    # Obtener valor de Pgd
    try:
        pgd_values = ampl.getVariable("Pgd").getValues().toPandas()
        pgd_values.reset_index(inplace=True)
        print("📊 Valores de la variable Pgd:")
        print(pgd_values)
        # Calcular HC
        HC = pgd_values["Pgd.val"].sum()  # Suma de los valores de Pgd
        print(f"✅ El valor de HC es: {HC}")
    except RuntimeError:
        print("⚠️ No se pudieron obtener valores para la variable 'Pgd'.")

    try:
        # Obtener los valores de V desde AMPL
        v_values = ampl.getVariable("V").getValues().toPandas()
        v_values.reset_index(inplace=True)
        #print("📊 Valores de la variable v:")
        #print(v_values)
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

    try:
        # Obtener los valores de PD y QD desde AMPL
        PD_values = ampl.getParameter("PD").getValues().toPandas()
        PD_values.reset_index(inplace=True)
        #print("📊 Valores del parametro PD:")
        #print(PD_values)
        QD_values = ampl.getParameter("QD").getValues().toPandas()
        QD_values.reset_index(inplace=True)
        #print("📊 Valores del parametro QD:")
        #print(QD_values)
    except RuntimeError:
            print("⚠️ No se pudieron obtener valores para los parametros 'PD' 'QD'.")

    return {
        "FuncionObjetivo": funcion_objetivo,
        "HC": HC,
        "V_min": v_min_pu,
        "V_max": v_max_pu,
        "Pgd_V": pgd_values,
        "PD_V": PD_values,
        "QD_V": QD_values
    }

if __name__ == "__main__":
    resultados = solve_optimization()
    #print("Resultados del problema de optimización:", resultados)


# Validación de resultados en las horas restantes
# Configuración de OpenDSS
dss_engine = dss.DSS
DSSText = dss_engine.Text
DSSCircuit = dss_engine.ActiveCircuit
DSSSolution = DSSCircuit.Solution
ControlQueue = DSSCircuit.CtrlQueue
dss_engine.AllowForms = 0

# Horas restantes que no estan dentro de las peores horas
horas_restantes = [hora for hora in range(1, 25) if hora not in peores_horas]
#print(horas_restantes)
#print(resultados.get("Pgd_V", None))
#print(resultados.get("PD_V", None))
#print(resultados.get("QD_V", None))

def generar_archivos_dss(hora, pgd_value, PD_value, QD_value):

    # Cargar perfiles de generación y demanda
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
    # 🔍 Filtrar por la hora específica
    generacion_hora = generacion_pivot.loc[:, hora]
    demanda_hora = demanda_pivot.loc[:, hora]

    # Combinar PD_value y QD_value en un solo DataFrame
    demanda_combined = PD_value.set_index('index').join(QD_value.set_index('index'), lsuffix='_PD', rsuffix='_QD')

    # 🔹 Generar archivo de generadores
    with open("IEEE_33_GENERADORES.dss", "w") as f:
        for _, row in pgd_value.iterrows():
            bus_number = int(row['index'])  
            perfil_generacion = generacion_hora[generacion_hora.index == bus_number].values
            pgd_value = row['Pgd.val'] * perfil_generacion[0] * 1000 if len(perfil_generacion) > 0 else 0
            f.write(f"New Generator.{bus_number} Bus1={bus_number} Phases=3 kV=12.66 kW={pgd_value:.6f} kvar=0\n")
    print(f"\n✅ Archivo de generadores para la hora {hora} generado con éxito.")

    # 🔹 Generar archivo de demandas
    with open("IEEE_33_DEMANDAS.dss", "w") as f:
        for _, row in demanda_combined.iterrows():
            bus_number = int(row.name)
            perfil_demanda = demanda_hora[demanda_hora.index == bus_number].values[0]
            PD_value = row['PD'] * perfil_demanda * 1000
            QD_value = row['QD'] * perfil_demanda * 1000
            f.write(f"New Load.{bus_number} Bus1={bus_number} Phases=3 kV=12.66 kW={PD_value:.6f} kvar={QD_value:.6f}\n")
    print(f"✅ Archivo de demandas para la hora {hora} generado con éxito.")


def obtener_resultados(horas):
    # Lista para almacenar las horas con tensiones fuera del rango
    horas_fuera_rango = []
    # Lista para almacenar las horas con tensiones fuera del rango y generación mayor a 0
    horas_fuera_rango_con_generacion = []

    # Extraer valores de parámetros en AMPL
    pgd_values = resultados["Pgd_V"]
    PD_values = resultados["PD_V"]
    QD_values = resultados["QD_V"]
    
    for hora in horas:
        generar_archivos_dss(hora, pgd_values, PD_values, QD_values)

        print(f"🔹 Ejecutando simulación en OpenDSS para la hora {hora}...")

        DSSText.Command = "Clear"
        DSSText.Command = "New Circuit.33bussystem bus1=1 basekv=12.66 pu=1 phases=3 mvasc3=1000000 frequency=60"
        if not DSSCircuit.Name:
            print(f"⚠️ Error: No se pudo crear el circuito en OpenDSS para la hora {hora}.")
        else:
            print(f"✅ Circuito creado correctamente para la hora {hora}.")
        with open("IEEE_33_GENERADORES.dss", "r") as file:
            comando_generadores = file.readlines()
        with open("IEEE_33_DEMANDAS.dss", "r") as file:
            comando_demandas = file.readlines()
        try:
            with open("IEEE_33_LINEAS.dss", "r") as file:
                comando_lineas = file.readlines()
        except FileNotFoundError:
            print("⚠️ El archivo 'IEEE_33_LINEAS.dss' no se encontró.")
            comando_lineas = []
        for generator in comando_generadores:
            generator = generator.strip()
            if generator:
                DSSText.Command = generator
        print(f"✅ {DSSCircuit.Generators.Count} Generadores cargados correctamente en OpenDSS para la hora {hora}.")
        for linea in comando_lineas:
            linea = linea.strip()
            if linea:
                DSSText.Command = linea
        print(f"✅ {DSSCircuit.Lines.Count} Líneas cargadas correctamente en OpenDSS para la hora {hora}.")
        for load in comando_demandas:
            load = load.strip()
            if load:
                DSSText.Command = load
        print(f"✅ {DSSCircuit.Loads.Count} Cargas cargadas correctamente en OpenDSS para la hora {hora}.")
        DSSText.Command = "Reset"
        DSSText.Command = "Set Voltagebases=[12.66]"
        DSSText.Command = "CalcVoltageBases"
        DSSText.Command = "Set ControlMode=Static"
        DSSText.Command = "Set Mode=Snapshot"
        DSSText.Command = "Set MaxControlIter=1000"
        DSSText.Command = "Solve"
        DSSSolution.Solve()
        if DSSSolution.Converged:
            print(f"✅ El circuito se resolvió correctamente para la hora {hora}")
        else:
            print(f"⚠️ La solución no converge para la hora {hora}")
        line_loss = DSSCircuit.LineLosses
        print("📢 Resultados pérdidas de potencia:")
        print(f"🔌Pérdidas Potencia Activa (P) = {round(abs(line_loss[0]), 3)} kW")
        print(f"🔌Pérdidas Potencia reactiva (Q) = {round(abs(line_loss[1]), 3)} kVAr")
        bus_names = []
        voltages_pu = []
        for bus in DSSCircuit.Buses:
            bus_names.append(bus.Name)
            voltages_pu.append(bus.puVmagAngle[0])
        max_voltage = max(voltages_pu)
        min_voltage = min(voltages_pu)
        # Obtener los nodos donde ocurren estos valores
        max_bus = bus_names[voltages_pu.index(max_voltage)]
        min_bus = bus_names[voltages_pu.index(min_voltage)]
        # Imprimir resultados
        print(f"⚡ Voltaje Máximo: {max_voltage:.4f} p.u. en el nodo {max_bus}")
        print(f"⚡ Voltaje Mínimo: {min_voltage:.4f} p.u. en el nodo {min_bus}")
        suma_generadores = DSSCircuit.Generators.kW
        print(suma_generadores)
        # Verificar si las tensiones están fuera del rango
        if max_voltage > 1.05 or min_voltage < 0.95:
            horas_fuera_rango.append(hora)
            # Verificar si la suma de los generadores instalados es mayor a 0
            if suma_generadores > 0:
                horas_fuera_rango_con_generacion.append(hora)
                print(f"⚠️ La hora {hora} tiene tensiones fuera del rango y es una hora de generación.")
            else:
                print(f"⚠️ La hora {hora} tiene tensiones fuera del rango pero no es una hora de generación.")

    # Imprimir las horas con tensiones fuera del rango
    print(f"\n⏰ Horas con tensiones fuera del rango 1.05 - 0.95: {horas_fuera_rango}")
    # Imprimir las horas con tensiones fuera del rango y generación mayor a 0
    print(f"⏰ Horas con tensiones fuera del rango y generación mayor a 0: {horas_fuera_rango_con_generacion}")
   # Unir las listas sin duplicados
    peores_horas_actualizadas = list(set(peores_horas + horas_fuera_rango_con_generacion))
    # Imprimir las horas combinadas
    print(f"⏰ Horas combinadas (peores horas y horas con generación): {peores_horas_actualizadas}")
    # Actualizar el archivo de peores horas
    horas_str = " ".join(map(str, peores_horas_actualizadas))  # Convertir horas a string para escribir en el archivo
    with open('IEEE_AMPL_PEORES_HORAS.dat', 'w') as f:
        f.write(f"set PH := {horas_str};\n")
    print(f"✅ Archivo de peores horas actualizado con las horas combinadas: {horas_str}")
# Ejecutar la validación para las horas restantes
obtener_resultados(horas_restantes)

