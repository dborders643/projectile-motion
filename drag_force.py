import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# ----------Setup---------- 
v0 = 10.0                           # launch speed [m/s]
phi_deg = 45.0                      # launch angle [deg]
phi = np.deg2rad(phi_deg)           # launch angle [rad]
vx0 = v0*np.cos(phi)                # launch speed (horizontal component) [m/s] 
vy0 = v0*np.sin(phi)                # launch speed (vertical component) [m/s]
t0 = 0.0                            # initial launch time [s]
x0 = 0.0                            # launch position (horizontal component) [m]
y0 = 0.0                            # launch position (vertical component) [m]
g = 9.81                            # gravational constant [m/s^2]
b = 0.1                             # drag coefficient [units don't matter]
tf = 3 * (2 * vy0/g)                # final time [s] (overshoot on purpose)
t = np.linspace(t0, tf, num=100)    # time array [s]
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

# Initial Conditions
x[0] = x0
y[0] = y0
vx[0] = vx0
vy[0] = vy0

# ----------Logic----------
for i in range(len(t) - 1):
    v = sqrt(vx**2 = vy**2)
    # Derivatives
    dxdt[i] = vx[i]
    dydt[i] = vy[i]
    dvxdt[i] = -b*v[i]**2*np.cos(phi)
    dvydt[i] = 
    # Ground Condition
    if y <= 0 and i > 0:
        break
    # Update
    dxdt[i+1] = dxdt[i] + dvxdt[i]*dt
    dydt[i+1] = dydt[i] + dvydt[i]*dt
    x[i+1] = x[i] + dxdt[i]*dt
    y[i+1] = y[i] + dydt[i]*dt
    t[i+1] = t[i] + dt
# ----------Visual---------- 
# Motion graph
plt.plot(x, y)
plt.xlabel("Horizontal Distance (m)")
plt.ylabel("Vertical Distance (m)")
plt.title("Projectile Motion with Drag")
plt.show()

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

ani = FuncAnimation(fig, update, frames=len(t),
                    init_func=init, interval=20)
plt.show()