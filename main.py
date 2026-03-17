import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

fig, ax = plt.subplots()

x = np.linspace(0, 2*np.pi, 100)
line, = ax.plot(x, np.sin(x))

def update(frame):
    line.set_ydata(np.sin(x + frame))  # update the data
    return line,

ani = FuncAnimation(
    fig,          # figure object
    update,       # function that updates each frame
    frames=100,   # number of frames
    interval=50   # milliseconds between frames
)

plt.show()