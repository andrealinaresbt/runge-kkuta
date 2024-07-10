import numpy as np
import matplotlib.pyplot as plt

# Parámetros del circuito
R = 29.08
L = 0.05
C = 10e-6
V_m = 110
omega = 377

# Ecuaciones diferenciales
def f(t, i, di):
    return (V_m * np.sin(omega * t) - R * i - (1/C) * i) / L

# Método de Runge-Kutta de 4º orden
def runge_kutta(t0, i0, di0, h, n):
    t = np.zeros(n)
    i = np.zeros(n)
    di = np.zeros(n)
    t[0] = t0
    i[0] = i0
    di[0] = di0
    for k in range(1, n):
        k1 = h * di[k-1]
        l1 = h * f(t[k-1], i[k-1], di[k-1])
        
        k2 = h * (di[k-1] + 0.5 * l1)
        l2 = h * f(t[k-1] + 0.5 * h, i[k-1] + 0.5 * k1, di[k-1] + 0.5 * l1)
        
        k3 = h * (di[k-1] + 0.5 * l2)
        l3 = h * f(t[k-1] + 0.5 * h, i[k-1] + 0.5 * k2, di[k-1] + 0.5 * l2)
        
        k4 = h * (di[k-1] + l3)
        l4 = h * f(t[k-1] + h, i[k-1] + k3, di[k-1] + l3)
        
        t[k] = t[k-1] + h
        i[k] = i[k-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
        di[k] = di[k-1] + (l1 + 2*l2 + 2*l3 + l4) / 6
        
    return t, i

# Parámetros de la simulación
t0 = 0
i0 = 0
di0 = 0
h = 1e-4  # Paso de tiempo más grande para mejorar eficiencia
n = int(0.02 / h)  # Simulación durante 20 ms

# Simulación
t, i = runge_kutta(t0, i0, di0, h, n)

# Cálculo del voltaje en R5
V_R5 = i * 10  # V_R5 = i * R5

# Cálculo del voltaje en L1
Z_L1 = 18.85  # Impedancia de L1
V_L1 = i * Z_L1  # V_L1 = i * Z_L1

# Ángulo de desfase del voltaje en L1
I_angle = -83.21  # Ángulo de la corriente calculado previamente
V_L1_angle = I_angle + 90  # Desfase total del voltaje en L1

# Gráficas en un mismo plot
plt.figure(figsize=(10, 6))
plt.plot(t, i, label='Corriente $i(t)$')
plt.plot(t, V_R5, label='Voltaje en $R5$ $V_{R5}(t)$')
plt.plot(t, V_L1, label='Voltaje en $L1$ $V_{L1}(t)$')

plt.title("Análisis del Circuito RLC")
plt.xlabel("Tiempo (s)")
plt.ylabel("Amplitud")
plt.legend()
plt.grid(True)
plt.show()

# Mostrar el ángulo de desfase del voltaje en L1
print(f"Ángulo de Desfase del Voltaje en L1: {V_L1_angle} grados")
