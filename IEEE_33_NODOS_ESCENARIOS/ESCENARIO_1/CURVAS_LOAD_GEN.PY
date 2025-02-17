## Curvas de Generación y Demanda
import matplotlib.pyplot as plt

# Curvas Demanda Diaria Sector Residencial E1 y E2
Hora = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
Demanda_Promedio_1 = [0.29, 0.31, 0.26, 0.27, 0.31, 0.33, 0.38, 0.37, 0.41, 0.47, 0.53, 0.55, 0.54, 0.50, 0.47, 0.42, 0.52, 0.61, 0.92, 1.00, 0.90, 0.75, 0.63, 0.45] # Estrato 1 y 2
Demanda_Promedio_2 = [0.30, 0.27, 0.26, 0.27, 0.24, 0.41, 0.44, 0.36, 0.36, 0.45, 0.41, 0.45, 0.56, 0.39, 0.41, 0.45, 0.47, 0.52, 0.77, 0.91, 0.98, 1.00, 0.64, 0.37] # Estrato 3
Demanda_Promedio_3 = [0.28, 0.25, 0.26, 0.25, 0.23, 0.44, 0.69, 0.44, 0.44, 0.36, 0.38, 0.51, 0.40, 0.38, 0.37, 0.41, 0.49, 0.34, 0.61, 0.81, 1.00, 0.86, 0.69, 0.36] # Estrato 4
Demanda_Promedio_4 = [0.45, 0.41, 0.38, 0.44, 0.39, 0.52, 0.68, 0.68, 0.59, 0.59, 0.63, 0.54, 0.59, 0.62, 0.60, 0.61, 0.48, 0.54, 0.72, 0.82, 1.00, 0.78, 0.67, 0.58] # Estrato 5
Demanda_Promedio_5 = [0.33, 0.38, 0.38, 0.34, 0.34, 0.36, 0.37, 0.57, 0.62, 0.52, 0.56, 1.00, 0.95, 0.74, 0.66, 0.65, 0.54, 0.48, 0.62, 0.73, 0.83, 0.73, 0.67, 0.55] # Estrato 6
Demanda_Promedio_6 = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.34, 0.76, 0.91, 0.93, 0.72, 0.41, 0.70, 0.95, 1.00, 0.98, 0.85, 0.61, 0.16, 0.07, 0.06, 0.05] # Comercial

plt.figure(figsize=(10, 6))
plt.plot(Hora, Demanda_Promedio_1, marker='o', label='Demanda 1')
plt.plot(Hora, Demanda_Promedio_2, marker='s', label='Demanda 2')
plt.plot(Hora, Demanda_Promedio_3, marker='^', label='Demanda 3')
plt.plot(Hora, Demanda_Promedio_4, marker='d', label='Demanda 4')
plt.plot(Hora, Demanda_Promedio_5, marker='x', label='Demanda 5')
plt.plot(Hora, Demanda_Promedio_6, marker='*', label='Demanda 5')
plt.xlabel("Hora del día")
plt.ylabel("Demanda Promedio")
plt.title("Curvas de Demanda Diaria - Sectores Residenciales")
plt.xticks(Hora)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()


# Curvas típicas de generación solar (normalizadas)
Generacion_Tipica_1 = [0, 0, 0, 0, 0, 0.10, 0.30, 0.50, 0.70, 0.85, 0.95, 1.00, 0.98, 0.92, 0.85, 0.70, 0.50, 0.30, 0.10,  0, 0, 0, 0, 0]  # Día despejado
Generacion_Tipica_2 = [0, 0, 0, 0, 0, 0.05, 0.20, 0.40, 0.60, 0.75, 0.85, 0.88, 0.85, 0.75, 0.60, 0.40, 0.20, 0.05,    0,  0, 0, 0, 0, 0]  # Día nublado
Generacion_Tipica_3 = [0, 0, 0, 0, 0, 0.08, 0.25, 0.45, 0.65, 0.82, 0.90, 0.94, 0.92, 0.88, 0.78, 0.62, 0.45, 0.20, 0.08,  0, 0, 0, 0, 0]  # Día con nubosidad parcial

plt.figure(figsize=(10, 6))
plt.plot(Hora, Generacion_Tipica_1, marker='o', label='Día despejado')
plt.plot(Hora, Generacion_Tipica_2, marker='s', label='Día nublado')
plt.plot(Hora, Generacion_Tipica_3, marker='^', label='Día con nubosidad parcial')
plt.xlabel("Hora del día")
plt.ylabel("Generación Solar Normalizada")
plt.title("Curvas Típicas de Generación Solar en Santander, Colombia")
plt.xticks(Hora)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()



