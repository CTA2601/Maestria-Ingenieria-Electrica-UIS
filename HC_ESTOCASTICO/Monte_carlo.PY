import dss
import os
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd

# --------------------- Parámetros configurables ---------------------
capacidad_max_gd = 1  # Capacidad máxima de los generadores en MW
num_iteraciones = 1000  # Número de iteraciones por nivel de penetración
paso_nivel_penetracion = 1  # Incremento del nivel de penetración en %
limite_flujo = 5000 # MW Limite de Flujo Inverso
# --------------------- Configuración de OpenDSS ---------------------
os.system('cls' if os.name == 'nt' else 'clear')
main_folder = "HC_ESTOCASTICO"
base_path = os.getcwd()
target_path = os.path.join(base_path, main_folder)
os.chdir(target_path)

print("📂 Directorio actual:", os.getcwd())

dss_engine = dss.DSS
if dss_engine.Start(0):
    print("✅ OpenDSS se inició correctamente.")
else:
    print("⚠️ Error al iniciar OpenDSS.")
    exit()

dss_text = dss_engine.Text
dss_circuit = dss_engine.ActiveCircuit
dss_solution = dss_circuit.Solution

# --------------------- Función para Generar Generadores ---------------------
def generar_archivo_generadores(nivel_penetracion):
    num_nodos = 33  
    num_gds = min(int((nivel_penetracion / 100) * num_nodos), num_nodos - 1)
    nodos_seleccionados = random.sample(range(2, num_nodos + 1), num_gds)
    generadores_config = {}

    with open("IEEE_33_GENERADORES.dss", "w") as f:
        for nodo in range(1, num_nodos + 1):  
            if nodo in nodos_seleccionados:
                capacidad = random.uniform(0, capacidad_max_gd * 1000)  
            else:
                capacidad = 0  
            
            config = f"New Generator.GD{nodo} Bus1={nodo} Phases=3 kV=12.66 kW={capacidad:.3f} kvar=0 Vminpu=0.7 Vmaxpu=1.2 conn=Wye Model=1\n"
            f.write(config)
            generadores_config[nodo] = config  

    return generadores_config  

# --------------------- Simulación Estocástica ---------------------
resultados = []
configuraciones_generadores = {}  
niveles_penetracion = list(range(0, 101, paso_nivel_penetracion))  

for nivel in niveles_penetracion:
    print(f"🔹 Nivel de penetración: {nivel}%")
    
    for _ in range(num_iteraciones):
        generadores_actual = generar_archivo_generadores(nivel)

        dss_text.Command = "Clear"
        dss_text.Command = "New Circuit.33bussystem bus1=1 basekv=12.66 pu=1 phases=3 mvasc3=20000000000 frequency=60"
        dss_text.Command = "Redirect IEEE_33_LINEAS.dss"
        dss_text.Command = "Redirect IEEE_33_DEMANDAS.dss"
        dss_text.Command = "Redirect IEEE_33_GENERADORES.dss"
        dss_text.Command = "New EnergyMeter.FeederMeter element=line.1-2 terminal=1"
        dss_text.Command = "Reset"
        dss_text.Command = "Set Voltagebases=[12.66]"
        dss_text.Command = "CalcVoltageBases"
        dss_text.Command = "Set ControlMode=OFF"
        dss_text.Command = "Set MaxControlIter=1000"
        dss_text.Command = "Solve" 
        dss_solution.Solve()

        if not dss_solution.Converged:
            continue

        line_loss = dss_circuit.LineLosses
        perdidas_p = round(abs(line_loss[0]), 3)  
        
        # Obtener potencia activa total desde OpenDSS
        total_active_power = dss_circuit.TotalPower[0]  # [0] es la potencia activa (MW)

        generator_names = dss_circuit.Generators.AllNames
        generator_power = []
        for gen in generator_names:
            dss_circuit.Generators.Name = gen
            generator_power.append(dss_circuit.Generators.kW)
        HC = sum(generator_power)

        # Obtener perfil de voltajes
        voltajes_pu = [bus.puVmagAngle[0] for bus in dss_circuit.Buses]
        max_voltage = max(voltajes_pu)
        min_voltage = min(voltajes_pu)

        # Verificar restricciones
        cumple_tension = 0.90 <= min_voltage <= 1.05 and 0.90 <= max_voltage <= 1.05
        flujo_subestacion = total_active_power  
        cumple_flujo = flujo_subestacion <= limite_flujo   

        cumple_limites = cumple_tension and cumple_flujo
        resultados.append((nivel, HC, perdidas_p, cumple_limites))

        if cumple_limites and (nivel not in configuraciones_generadores or HC > configuraciones_generadores[nivel][0]):  
            configuraciones_generadores[nivel] = (HC, generadores_actual)

# --------------------- Filtrar Resultados ---------------------
df_resultados = pd.DataFrame(resultados, columns=["Nivel", "HC", "Perdidas", "Cumple"])
validos = df_resultados[df_resultados["Cumple"] == True]

# Obtener valores clave
max_HC = validos["HC"].max()
primer_invalido = df_resultados[df_resultados["Cumple"] == False]
HC_primer_invalido = primer_invalido["HC"].iloc[0] if not primer_invalido.empty else None

# Obtener valores clave
max_HC_nivel = validos.loc[validos["HC"].idxmax(), "Nivel"]
nivel_min_perdida = validos.loc[validos["Perdidas"].idxmin(), "Nivel"]

# 🔹 Guardar HC_MAX.dss con generadores en todos los nodos
with open("HC_MAX.dss", "w") as f:
    for nodo in range(1, 34):  
        f.write(configuraciones_generadores[max_HC_nivel][1].get(nodo, f"New Generator.GD{nodo} Bus1={nodo} Phases=3 kV=12.66 kW=0 kvar=0 Vminpu=0.7 Vmaxpu=1.2 conn=Wye Model=1\n"))

# 🔹 Guardar MIN_LOSE.dss con generadores en todos los nodos
with open("MIN_LOSE.dss", "w") as f:
    for nodo in range(1, 34):
        f.write(configuraciones_generadores[nivel_min_perdida][1].get(nodo, f"New Generator.GD{nodo} Bus1={nodo} Phases=3 kV=12.66 kW=0 kvar=0 Vminpu=0.7 Vmaxpu=1.2 conn=Wye Model=1\n"))


# --------------------- Gráfica de HC vs Nivel de Penetración ---------------------
plt.figure(figsize=(12, 6))
plt.scatter(df_resultados["Nivel"], df_resultados["HC"], 
            c=df_resultados["Cumple"].map({True: "green", False: "red"}), 
            alpha=0.6, edgecolors='k')

# Líneas de referencia
plt.axhline(y=max_HC, color='blue', linestyle="--", label=f"Máximo HC Válido: {max_HC:.2f} kW")
if HC_primer_invalido is not None:
    plt.axhline(y=HC_primer_invalido, color='purple', linestyle="--", label=f"Primer HC Inválido: {HC_primer_invalido:.2f} kW")

plt.xlabel("Nivel de penetración (%)")
plt.ylabel("Hosting Capacity (HC) [kW]")
plt.title("Hosting Capacity vs Nivel de Penetración")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.show()

# --------------------- Gráfica de Pérdidas vs HC ---------------------
plt.figure(figsize=(12, 6))
plt.scatter(validos["HC"], validos["Perdidas"], c="blue", alpha=0.6, edgecolors='k')

min_perdida = validos["Perdidas"].min()
hc_min_perdida = validos.loc[validos["Perdidas"].idxmin(), "HC"]
plt.axhline(y=min_perdida, color='blue', linestyle="--", label=f"Mínimas pérdidas: {min_perdida:.2f} kW en HC={hc_min_perdida:.2f} kW")

plt.xlabel("Hosting Capacity (HC) [kW]")
plt.ylabel("Pérdidas de energía (kW)")
plt.title("Pérdidas de Energía vs Hosting Capacity (Solo HC Válidos)")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.show()
