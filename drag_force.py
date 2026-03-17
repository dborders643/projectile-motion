import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# ----------Setup---------- 
v0 = 100.0                          # launch speed [m/s]
phi_deg = 45.0                      # launch angle [deg]
phi = np.deg2rad(phi_deg)           # launch angle [rad]
vx0 = v0*np.cos(phi)                # launch speed (horizontal component) [m/s] 
vy0 = v0*np.sin(phi)                # launch speed (vertical component) [m/s]
t0 = 0.0                            # initial launch time [s]
x0 = 0.0                            # launch position (horizontal component) [m]
y0 = 0.0                            # launch position (vertical component) [m]
g = 9.81                            # gravational constant [m/s^2]
b = 0.1                             # drag coefficient [units don't matter]
tf = 100                            # final time [s] (overshoot on purpose)
t = np.linspace(t0, tf, num=10000)  # time array [s]
dt = t[1] - t[0]                    # scalar time step [s]

# Preallocation
x = np.zeros_like(t)
y = np.zeros_like(t)
dxdt = np.zeros_like(t)
dydt = np.zeros_like(t)
dvxdt = np.zeros_like(t)
dvydt = np.zeros_like(t)
vx = np.zeros_like(t)
vy = np.zeros_like(t)
v = np.zeros_like(t)

# Initial Conditions
x[0] = x0
y[0] = y0
vx[0] = vx0
vy[0] = vy0
v[0] = v0

# ----------Logic----------
for i in range(len(t) - 1):
    # Derivatives
    dxdt[i] = vx[i]
    dydt[i] = vy[i]
    v[i] = np.sqrt(vx[i]**2 + vy[i]**2)
    dvxdt[i] = -b*v[i]*vx[i]
    dvydt[i] = -g - b*v[i]*vy[i]
    # Ground Condition
    if y[i] <= 0 and i > 0:
        break
    # Update
    vx[i+1] = vx[i] + dvxdt[i]*dt
    vy[i+1] = vy[i] + dvydt[i]*dt
    x[i+1] = x[i] + dxdt[i]*dt
    y[i+1] = y[i] + dydt[i]*dt
    t[i+1] = t[i] + dt
# Trim zeros
x = x[:i]
y = y[:i]

# ----------Animation----------
fig, ax = plt.subplots()
ax.plot(x, y, 'k--', alpha=0.3)
ball, = ax.plot([], [], 'ro', markersize=8)

def init():
    ball.set_data([], [])
    return ball,

def update(i):
    ball.set_data([x[i]], [y[i]])
    return ball,

ani = FuncAnimation(fig, update, frames=len(x),
                    init_func=init, interval=20)
ax.set_xlabel("Horizontal Distance (m)")
ax.set_ylabel("Vertical Distance (m)")
ax.set_title("Projectile Motion with Drag")
plt.show()