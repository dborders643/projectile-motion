import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# ----------Setup----------
g = 9.81                            # gravitational constant [m/s^2]
l = 1.0                             # pendulum length [m]
q = 0.5                             # damping coefficient [1/s]
wD = 1.0                            # driving force frequency [rad/s]
FD = 1.2                            # driving force amplitude [rad/s^2]
phi0 = np.pi/4                      # initial angle [rad]
w0 = 0.0                            # initial angular velocity [rad/s]
t0 = 0.0                            # initial time [s]
tf = 10.0                           # final time [s]
t = np.linspace(t0, tf, num=1501)   # time array [s]
dt = t[1] - t[0]                    # scalar time step [s]

# Preallocation
phi = np.zeros_like(t)
w = np.zeros_like(t)
dwdt = np.zeros_like(t)

# ----------Initial Conditions----------
phi[0] = phi0
w[0] = w0

# ----------Logic----------
for i in range(len(t) - 1):
    # Derivatives 
    dwdt[i] = -(g/l)*np.sin(phi[i]) - q*w[i] + FD*np.sin(wD*t[i])
    # Update
    w[i+1] = w[i] + dwdt[i]*dt
    phi[i+1] = phi[i] + w[i+1]*dt
    # Angular Position Condition
    if phi[i+1] < -np.pi:
        phi[i+1] += 2*np.pi
    elif phi[i+1] > np.pi:
        phi[i+1] -= 2*np.pi

# ----------Animation----------
# Convert to cartesian coordinates
px = l * np.sin(phi)
py = -l * np.cos(phi)

# Layout
fig = plt.figure(figsize=(12, 6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.5, 1])
ax_anim = fig.add_subplot(gs[:, 0])       # left, full height
ax_time = fig.add_subplot(gs[0, 1])       # top right
ax_phase = fig.add_subplot(gs[1, 1])      # bottom right

# Animation axes
ax_anim.set_xlim(-l*1.2, l*1.2)
ax_anim.set_ylim(-l*1.2, l*1.2)
ax_anim.set_aspect('equal')
ax_anim.set_title('Pendulum')
ax_anim.set_xlabel('x (m)')
ax_anim.set_ylabel('y (m)')
rod, = ax_anim.plot([], [], 'b-', linewidth=2)
bob, = ax_anim.plot([], [], 'ro', markersize=10)
time_text = ax_anim.text(0.05, 0.95, '', transform=ax_anim.transAxes)

# Time series — empty initial lines
line_phi, = ax_time.plot([], [], 'r-', label='θ (rad)')
line_w, = ax_time.plot([], [], 'b-', label='ω (rad/s)')
ax_time.set_xlim(t0, tf)
ax_time.set_ylim(min(min(phi), min(w))*1.2, max(max(phi), max(w))*1.2)
ax_time.set_xlabel('Time (s)')
ax_time.set_ylabel('θ (rad), ω (rad/s)')
ax_time.set_title('Angular Position and Velocity')
ax_time.legend()

# Phase space — empty initial line
line_phase, = ax_phase.plot([], [], 'g-', linewidth=0.5)
ax_phase.set_xlim(-np.pi, np.pi)
ax_phase.set_ylim(min(w)*1.2, max(w)*1.2)
ax_phase.set_xlabel('θ (rad)')
ax_phase.set_ylabel('ω (rad/s)')
ax_phase.set_title('Phase Space')

def init():
    rod.set_data([], [])
    bob.set_data([], [])
    time_text.set_text('')
    line_phi.set_data([], [])
    line_w.set_data([], [])
    line_phase.set_data([], [])
    return rod, bob, time_text, line_phi, line_w, line_phase,

def update(i):
    rod.set_data([0, px[i]], [0, py[i]])
    bob.set_data([px[i]], [py[i]])
    time_text.set_text(f't = {t[i]:.2f} s')
    line_phi.set_data(t[:i], phi[:i])
    line_w.set_data(t[:i], w[:i])
    line_phase.set_data(phi[:i], w[:i])
    return rod, bob, time_text, line_phi, line_w, line_phase,

ani = FuncAnimation(fig, update, frames=range(0, len(t), 10),
                    init_func=init, interval=1)
plt.show()