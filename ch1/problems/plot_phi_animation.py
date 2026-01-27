import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = np.loadtxt('phi_history.csv', skiprows=1, delimiter=',')
iters = data[:,0].astype(int)
phi_vals = data[:,1:]

n_nodes = phi_vals.shape[1]
x = np.linspace(0, 0.1, n_nodes)

fig, ax = plt.subplots(figsize=(8,5))
line, = ax.plot(x, phi_vals[0], '-o', markersize=4, color='k')
ax.set_xlim(x[0], x[-1])
ax.set_xlabel('x (m)')
ax.set_ylabel('phi (V)')
title = ax.set_title(f'iter 0')
ax.grid(True)

def update(frame):
    y = phi_vals[frame]
    line.set_ydata(y)
    ymin, ymax = y.min(), y.max()
    margin = (ymax - ymin) * 0.1 if ymax != ymin else 1
    ax.set_ylim(ymin - margin, ymax + margin)
    title.set_text(f'iter {iters[frame]}')
    return line, title

ani = FuncAnimation(fig, update, frames=len(iters), interval=500, blit=False)
plt.tight_layout()
plt.show()
