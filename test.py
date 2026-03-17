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
y = np.zeros_like(dt)

# Equations
for _ in range(np.size(dt)):
    
# Visualize 
