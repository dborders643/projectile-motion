import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# ----------Setup----------
v = 1.0                             # wave speed [m/s]
L = 1.0                             # string length [m]
t0 = 0.0                            # initial time [s]
tf = 10.0                           # final time [s]
nx = 200                            # number of spatial points
nt = 1000                           # number of time points
x = np.linspace(0, L, nx)          # spatial array [m]
t = np.linspace(t0, tf, nt)        # time array [s]
dx = x[1] - x[0]                   # spatial step [m]
dt = t[1] - t[0]                   # time step [s]

# Courant number -- must be <= 1 for stability!
C = v * dt / dx
print(f'Courant number: {C:.3f}')

# ----------Preallocation----------
# u is a 2D array -- rows are time steps, columns are spatial points
# TODO: preallocate u with shape (nt, nx) using np.zeros

# ----------Initial Conditions----------
# A gaussian pulse in the middle of the string
u[0, :] = np.exp(-((x - L/2)**2) / 0.01)  # initial shape
u[1, :] = u[0, :]                           # initial velocity = 0

# ----------Boundary Conditions----------
# Fixed ends (string pinned at both ends)
# u[:, 0] = 0 and u[:, -1] = 0 (already zero from np.zeros)

# ----------Logic----------
# TODO: fill in the finite difference loop
# loop over time from i=1 to nt-2
# for each i, loop over space from j=1 to nx-2 (skip boundaries!)
#   u[i+1, j] = 2*u[i,j] - u[i-1,j] + C**2 * (u[i,j+1] - 2*u[i,j] + u[i,j-1])

# ----------Animation----------
# TODO: animate u[i, :] vs x for each time step i