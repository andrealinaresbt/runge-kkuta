import numpy as np
import matplotlib.pyplot as plt
import math

# Parámetros del circuito
R = 29.08           # Resistencia en ohmios
L = 0.15            # Inductancia total en henrios (L1 + L2)
C = 10*10^-6           # Capacitancia en faradios
V_m = 55            # Amplitud de la fuente de voltaje en voltios
omega = 120 * math.pi  # Frecuencia angular en radianes por segundo

# Impedancia total del circuito
Z_total = R + 1j*omega*L - 1j/(omega*C)
print(f"Z_total = {Z_total.real:.2f} + {Z_total.imag:.2f}j Ω")

# Definimos las ecuaciones diferenciales del sistema
def f(t, x):
    i, di = x
    d2i_dt2 = (V_m * omega * np.cos(omega * t) - R * di - (1/C) * i) / L
    return np.array([di, d2i_dt2])

# Método de Runge-Kutta de 4º orden
def runge_kutta(t0, x0, h, n):
    t = np.zeros(n)
    x = np.zeros((n, len(x0)))
    
    t[0] = t0
    x[0] = x0
    
    for k in range(1, n):
        k1 = h * f(t[k-1], x[k-1])
        k2 = h * f(t[k-1] + 0.5*h, x[k-1] + 0.5*k1)
        k3 = h * f(t[k-1] + 0.5*h, x[k-1] + 0.5*k2)
        k4 = h * f(t[k-1] + h, x[k-1] + k3)
        
        t[k] = t[k-1] + h
        x[k] = x[k-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return t, x

# Parámetros de la simulación
t0 = 0                  # Tiempo inicial
i0 = 0                  # Corriente inicial
di0 = 0                 # Derivada de la corriente inicial
h = 1e-4                # Paso de tiempo (ajustado para eficiencia)
n = int(0.1 / h)        # Número de pasos (simulación durante 100 ms)

# Simulación
t, x = runge_kutta(t0, [i0, di0], h, n)
i = x[:, 0]

# Cálculo del voltaje en R5
R5 = 10
V_R5 = i * R5  # V_R5 = i * R5

# Cálculo del voltaje en L1
Z_L1 = 18.85  # Impedancia de L1
V_L1 = i * Z_L1  # V_L1 = i * Z_L1

# Ángulo de desfase del voltaje en L1
I_angle = np.angle(np.fft.fft(i)[1], deg=True)  # Ángulo de la corriente calculado previamente
V_L1_angle = I_angle + 90  # Desfase total del voltaje en L1

# Funciones de salida
i_func = lambda t: V_m / np.abs(Z_total) * np.sin(omega * t - np.angle(Z_total))
V_R5_func = lambda t: i_func(t) * R5
V_L1_func = lambda t: i_func(t) * Z_L1

# Mostrar las funciones
print(f"i(t) = {V_m / np.abs(Z_total):.2f} * sin({omega} * t - {np.angle(Z_total):.2f})")
print(f"V_R5(t) = {V_m / np.abs(Z_total) * R5:.2f} * sin({omega} * t - {np.angle(Z_total):.2f})")
print(f"V_L1(t) = {V_m / np.abs(Z_total) * Z_L1:.2f} * sin({omega} * t - {np.angle(Z_total):.2f})")

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
print(f"Ángulo de Desfase del Voltaje en L1: {V_L1_angle:.2f} grados")
