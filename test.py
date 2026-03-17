import matplotlib.pyplot as plt
import numpy as np

# Givens
v0 = 10.0       # launch speed [m/s]
phi = 45.0      # launch angle [deg]
t0 = 0.0        # initial launch time [s]
t = 10.0        # final time [s]
x0 = 0.0        # launch position (horizontal component) [m]
y0 = 0.0        # launch position (vertical component) [m]

# Setup 
dt = np.linspace(t0, t, num=1000)
x = np.zeros_like(dt)
# Preallocation
y = np.zeros_like(dt)
vx = np.zeros_like(dt)
vy = np.zeros_like(dt)

vx0 = v0*np.cos(phi)
vy0 = v0*np.sin(phi)

np.cos()

# Equations
for i in range(np.size(dt)):
    

    # Update Equations (Euler Method)
    x(i+1) = x(i) + dvxdt(i)*dt(i)
    y(i+1) = y(i) + dvydt(i)*dt(i)
# Visualize 
