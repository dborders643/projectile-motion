import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# ----------Setup----------
v = 1.0                             # wave speed [m/s]
L = 1.0                             # string length [m]
tf = 10.0                           # final time [s]
nx = 200                            # number of spatial points
nt = 2000                           # number of time points
x = np.linspace(0, L, nx)           # spatial array [m]
t = np.linspace(0, tf, nt)          # time array [s]
dx = x[1] - x[0]                    # spatial step [m]
dt = t[1] - t[0]                    # time step [s]

# Courant number -- must be <= 1 for stability
C = v * dt / dx
print(f'Courant number: {C:.3f}')

# Preallocation
u = np.zeros((nt, nx))

# Initial Conditions
# A gaussian pulse in the middle of the string
u[0, :] = np.exp(-((x - L/2)**2) / 0.01)
u[1, :] = u[0, :]  # v0 = 0

# ----------Logic----------
for i in range(1, nt-1):
    for j in range(1, nx-1):
        u[i+1, j] = 2*u[i, j] - u[i-1, j] + C**2*(u[i, j+1] - 2*u[i, j] + u[i, j-1])

# ----------Animation----------
fig, ax = plt.subplots()
line, = ax.plot(x, u[0, :], 'b-', linewidth=2)
ax.set_xlim(0, L)
ax.set_ylim(-1.2, 1.2)
ax.set_xlabel('Position (m)')
ax.set_ylabel('Displacement (m)')
ax.set_title('1D Wave Equation')
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)

def init():
    line.set_ydata(u[0, :])
    time_text.set_text('')
    return line, time_text,

def update(i):
    line.set_ydata(u[i, :])
    time_text.set_text(f't = {t[i]:.2f} s')
    return line, time_text,

ani = FuncAnimation(fig, update, frames=range(0, nt, 2),
                    init_func=init, interval=20)
plt.show()