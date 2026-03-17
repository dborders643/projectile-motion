import matplotlib.pyplot as plt
import numpy as np

# ----------Setup---------- 
v0 = 10.0                          # launch speed [m/s]
phi_deg = 45.0                      # launch angle [deg]
phi = np.deg2rad(phi_deg)           # launch angle [rad]
vx0 = v0*np.cos(phi)                # launch speed (horizontal component) [m/s] 
vy0 = v0*np.sin(phi)                # launch speed (vertical component) [m/s]
t0 = 0.0                            # initial launch time [s]
x0 = 0.0                            # launch position (horizontal component) [m]
y0 = 0.0                            # launch position (vertical component) [m]
g = 9.81                            # gravational constant [m/s^2]
tf = 2 * vy0/g                      # final time [s]
t = np.linspace(t0, tf, num=100)    # time array [s]

# ----------Logic----------
vx = vx0
vy = vy0 - g*t
x = x0 + vx0*t
y = y0 + vy0*t - 0.5*g*t**2
# ----------Debug----------
print(f'x values: {x[0:6]}')
print(f'y values: {y[0:6]}')
# ----------Visual---------- 
plt.plot(x, y)
plt.xlabel("Horizontal Distance (m)")
plt.ylabel("Vertical Distance (m)")
plt.title("Projectile Motion")
plt.show()