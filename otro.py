import numpy as np
import math

# Parámetros del circuito
R1 = 15
R2 = 10
R3 = 10
R4 = 5
R5 = 10
L1 = 0.05
L2 = 0.1
C = 10 *10^-6
omega = 120 * math.pi

# Calcular Z_{R2R3}
Z_R2R3 = (R2 * R3) / (R2 + R3)

# Calcular Z_{R4L2}
Z_R4L2 = R4 + 1j * omega * L2

# Calcular Z_{R5R4L2}
Z_R5R4L2 = (R5 * Z_R4L2) / (R5 + Z_R4L2)

# Calcular Z_{R2R3R5R4L2}
Z_R2R3R5R4L2 = Z_R2R3 + Z_R5R4L2

# Calcular Z_{C}
Z_C = -1j / (omega * C)

# Calcular Z_total
Z_total = R1 + 1j * omega * L1 + Z_R2R3R5R4L2 + Z_C

# Mostrar el resultado
print(f"Z_total = {Z_total.real:.2f} + {Z_total.imag:.2f}j Ω")