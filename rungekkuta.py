import numpy as np
import matplotlib.pyplot as plt

# Parámetros del circuito
R = 29.08           # Resistencia en ohmios
L = 0.05            # Inductancia en henrios
C = 10e-6           # Capacitancia en faradios
V_m = 110           # Amplitud de la fuente de voltaje en voltios
omega = 377         # Frecuencia angular en radianes por segundo

# Ecuaciones diferenciales
def f(t, i, di):
    """
    Función que define la ecuación diferencial del circuito RLC.
    Calcula la derivada de la corriente respecto al tiempo.

    :param t: Tiempo actual
    :param i: Corriente actual
    :param di: Derivada de la corriente actual
    :return: Derivada de la corriente respecto al tiempo
    """
    return (V_m * np.sin(omega * t) - R * i - (1/C) * i) / L

# Método de Runge-Kutta de 4º orden
def runge_kutta(t0, i0, di0, h, n):
    """
    Implementación del método de Runge-Kutta de 4º orden para resolver la ecuación diferencial del circuito RLC.

    :param t0: Tiempo inicial
    :param i0: Valor inicial de la corriente
    :param di0: Valor inicial de la derivada de la corriente
    :param h: Paso de tiempo
    :param n: Número de pasos
    :return: Arrays de tiempo y corriente calculados
    """
    t = np.zeros(n)   # Array para almacenar los valores de tiempo
    i = np.zeros(n)   # Array para almacenar los valores de corriente
    di = np.zeros(n)  # Array para almacenar los valores de la derivada de la corriente

    t[0] = t0        # Inicialización del tiempo
    i[0] = i0        # Inicialización de la corriente
    di[0] = di0      # Inicialización de la derivada de la corriente

    for k in range(1, n):
        # Calcular los coeficientes k y l (k1, l1, k2, l2, k3, l3, k4, l4)
        k1 = h * di[k-1]
        l1 = h * f(t[k-1], i[k-1], di[k-1])
        
        k2 = h * (di[k-1] + 0.5 * l1)
        l2 = h * f(t[k-1] + 0.5 * h, i[k-1] + 0.5 * k1, di[k-1] + 0.5 * l1)
        
        k3 = h * (di[k-1] + 0.5 * l2)
        l3 = h * f(t[k-1] + 0.5 * h, i[k-1] + 0.5 * k2, di[k-1] + 0.5 * l2)
        
        k4 = h * (di[k-1] + l3)
        l4 = h * f(t[k-1] + h, i[k-1] + k3, di[k-1] + l3)
        
        # Actualizar los valores de tiempo, corriente y derivada de la corriente
        t[k] = t[k-1] + h
        i[k] = i[k-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
        di[k] = di[k-1] + (l1 + 2*l2 + 2*l3 + l4) / 6
        
    return t, i

# Parámetros de la simulación
t0 = 0              # Tiempo inicial
i0 = 0              # Corriente inicial
di0 = 0             # Derivada de la corriente inicial
h = 1e-4            # Paso de tiempo (ajustado para eficiencia)
n = int(0.02 / h)   # Número de pasos (simulación durante 20 ms)

# Simulación
t, i = runge_kutta(t0, i0, di0, h, n)

# Cálculo del voltaje en R5
V_R5 = i * 10  # V_R5 = i * R5

# Cálculo del voltaje en L1
Z_L1 = 18.85  # Impedancia de L1
V_L1 = i * Z_L1  # V_L1 = i * Z_L1

# Ángulo de desfase del voltaje en L1
I_angle = 83.21  # Ángulo de la corriente calculado previamente
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
