import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# ----------Setup----------
v = 1.0                                 # wave speed [m/s]
Lx = 1.0                                # x length [m]
Ly = 1.0                                # y length [m]
tf = 10.0                               # final time [s]
nx = 100                                # number of x spatial points
ny = 100                                # number of y spatial points
nt = 4000                               # number of time points
x = np.linspace(0, Lx, nx)              # x spatial array [m]
y = np.linspace(0, Ly, ny)              # y spatial array [m]
t = np.linspace(0, tf, nt)              # time array [s]
dx = x[1] - x[0]                        # x spatial step [m]
dy = y[1] - y[0]                        # y spatial step [m]
dt = t[1] - t[0]                        # time step [s]

# Courant number -- must be <= 1 for stability
C = v * dt / dx
print(f'Courant number: {C:.3f}')

# ----------Preallocation----------
u = np.zeros((nt, nx, ny))

# ----------Initial Conditions----------
# Gaussian pulse in the center of the membrane
xx, yy = np.meshgrid(x, y, indexing='ij')
u[0, :, :] = np.exp(-((xx - Lx/2)**2 + (yy - Ly/2)**2) / 0.01)
u[1, :, :] = u[0, :, :] # v0 = 0

# ----------Logic----------
for i in range(1, nt - 1):
    u[i+1, 1:-1, 1:-1] = (2*u[i, 1:-1, 1:-1] - u[i-1, 1:-1, 1:-1]
        + C**2 * (u[i, 2:, 1:-1] - 2*u[i, 1:-1, 1:-1] + u[i, :-2, 1:-1])
        + C**2 * (u[i, 1:-1, 2:] - 2*u[i, 1:-1, 1:-1] + u[i, 1:-1, :-2]))

# ----------Animation----------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.set_zlim(-1.0, 1.0)
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('Displacement (m)')
ax.set_title('2D Wave Equation')

surf = [ax.plot_surface(xx, yy, u[0, :, :], cmap='coolwarm')]

def update(i):
    surf[0].remove()
    surf[0] = ax.plot_surface(xx, yy, u[i, :, :], cmap='coolwarm')
    ax.set_title(f't = {t[i]:.2f} s')
    return surf

ani = FuncAnimation(fig, update, frames=range(0, nt, 5), interval=1)
plt.show()