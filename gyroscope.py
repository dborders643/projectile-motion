import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import numpy as np

# ----------Setup----------
m       = 1.0                               # disc mass [kg]
R       = 0.5                               # disc radius [m]
d       = 0.6                               # pivot-to-COM distance [m]
omega   = 12.0                              # spin rate [rad/s]
g       = 9.81                              # gravitational constant [m/s^2]
theta   = np.deg2rad(25.0)                  # tilt of spin axis from vertical [rad]

I       = 0.5 * m * R**2                    # moment of inertia (solid disc) [kg m^2]
L_mag   = I * omega                         # magnitude of angular momentum [kg m^2/s]
tau_mag = m * g * d * np.sin(theta)         # magnitude of torque [N m]
Omega_p = tau_mag / L_mag                   # precession rate [rad/s]

t0      = 0.0                               # initial time [s]
tf      = 2 * np.pi / Omega_p               # one full precession period [s]
t       = np.linspace(t0, tf, num=500)      # time array [s]

# ----------Logic----------
phi     = Omega_p * t                       # precession angle over time [rad]

# Spin axis unit vector as a function of precession angle
spin_x  = np.sin(theta) * np.cos(phi)      # x-component of spin axis
spin_y  = np.sin(theta) * np.sin(phi)      # y-component of spin axis
spin_z  = np.cos(theta) * np.ones_like(t)  # z-component of spin axis

# COM position (pivot at origin, COM along spin axis at distance d)
com_x   = d * spin_x                        # COM x [m]
com_y   = d * spin_y                        # COM y [m]
com_z   = d * spin_z                        # COM z [m]

# Angular momentum vector L = I*omega * spin_axis
L_x     = L_mag * spin_x                   # L x-component [kg m^2/s]
L_y     = L_mag * spin_y                   # L y-component [kg m^2/s]
L_z     = L_mag * spin_z                   # L z-component [kg m^2/s]

# Torque vector tau = r x F,  F = [0, 0, -mg]
# r = [com_x, com_y, com_z], F = [0, 0, -mg]
# tau = (-mg*com_y, mg*com_x, 0)
tau_x   = -m * g * com_y                   # torque x-component [N m]
tau_y   =  m * g * com_x                   # torque y-component [N m]
tau_z   =  np.zeros_like(t)                # torque z-component [N m]

# Linear momentum of COM:  p = m * v_com,  v_com = Omega_p_vec x r_com
# Omega_p_vec = [0, 0, Omega_p]  =>  v = (-Omega_p*cy, Omega_p*cx, 0)
p_x     = m * (-Omega_p * com_y)           # p x-component [kg m/s]
p_y     = m * ( Omega_p * com_x)           # p y-component [kg m/s]
p_z     = np.zeros_like(t)                 # p z-component [kg m/s]

# Force (gravity, constant downward)
F_x     = np.zeros_like(t)                 # F x-component [N]
F_y     = np.zeros_like(t)                 # F y-component [N]
F_z     = -m * g * np.ones_like(t)         # F z-component [N]

# ----------Animation----------
fig     = plt.figure(figsize=(8, 7))
ax      = fig.add_subplot(111, projection='3d')

# Ghost precession cone trace (drawn once, like the dashed arc in projectile)
ax.plot(com_x, com_y, com_z, 'k--', alpha=0.2)

# Artists updated each frame
rod,    = ax.plot([], [], [], color='#888888', linewidth=2.5)
disc,   = ax.plot([], [], [], color='#4A9EFF', linewidth=1.5)

# Quivers cannot be updated in place, so track and remove each frame
quivers = []

SCALE   = 0.07     # vector arrow display scale

def init():
    rod.set_data_3d([], [], [])
    disc.set_data_3d([], [], [])
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-0.1, 1.4)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    ax.set_title('Gyroscopic Precession')
    return rod, disc

def update(i):
    for q in quivers:
        q.remove()
    quivers.clear()

    cx, cy, cz = com_x[i], com_y[i], com_z[i]

    # Pivot rod
    rod.set_data_3d([0, cx], [0, cy], [0, cz])

    # Disc rim (circle perpendicular to spin axis)
    n    = np.array([spin_x[i], spin_y[i], spin_z[i]])
    perp = np.array([1,0,0]) if abs(spin_x[i]) < 0.9 else np.array([0,1,0])
    u    = np.cross(n, perp); u /= np.linalg.norm(u)
    v    = np.cross(n, u)
    ang  = np.linspace(0, 2*np.pi, 64)
    disc.set_data_3d(
        cx + R*(np.cos(ang)*u[0] + np.sin(ang)*v[0]),
        cy + R*(np.cos(ang)*u[1] + np.sin(ang)*v[1]),
        cz + R*(np.cos(ang)*u[2] + np.sin(ang)*v[2])
    )

    # All five vectors: (origin, direction, scale, color, label)
    vectors = [
        ([cx,cy,cz], [L_x[i],   L_y[i],   L_z[i]  ], SCALE,      "#2F5A08", 'L'),
        ([cx,cy,cz], [tau_x[i], tau_y[i],  tau_z[i]], SCALE*3,    '#FF6B4A', 'tau'),
        ([0, 0, 0 ], [cx,       cy,        cz      ], 1.0,        '#2EC28A', 'r'),
        ([cx,cy,cz], [p_x[i],   p_y[i],    p_z[i]  ], SCALE*0.5,  '#F5A623', 'p'),
        ([cx,cy,cz], [F_x[i],   F_y[i],    F_z[i]  ], SCALE*0.5,  '#C06FE0', 'F'),
    ]
    for (ox,oy,oz), (dx,dy,dz), scale, color, label in vectors:
        q = ax.quiver(ox, oy, oz, dx*scale, dy*scale, dz*scale,
                      color=color, linewidth=2, arrow_length_ratio=0.2, label=label)
        quivers.append(q)

    ax.legend(loc='upper left', fontsize=8, framealpha=0.5)
    return rod, disc

ani = FuncAnimation(
    fig,
    update,
    frames=range(1, len(t), 5),
    init_func=init,
    interval=1)

plt.tight_layout()
plt.show()