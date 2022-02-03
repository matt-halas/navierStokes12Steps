# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython import display

n_points = 41
domainLength = 2
x_loc = np.linspace(0, domainLength, n_points)
dx = domainLength / n_points
dt = 0.025
n_time_steps = 25
c = 1

def addWave(u, start, end):
    idx_start = int(len(u) * start)
    idx_end = int(len(u) * end)
    u[idx_start:idx_end] += 1
    return u

# %%
def timeStep(u, c, dx, dt):
    u_diff = u[1:] - u[:-1]
    u_sub = u_diff[0]
    u_diff = np.append(u_sub, u_diff)
    u_new = u - c * dt / dx * u_diff
    return u_new

u = np.ones((n_points))
u = addWave(u, 0.25, 0.5)
u_all = np.zeros((n_time_steps, n_points))

for i in range(n_time_steps):
    u = timeStep(u, c, dx, dt)
    u_all[i, :] = u

def plotVelocity(i):
    if i >= len(u_all):
        plt.close()
    else:
        plt.clf()
        plt.xlim((0, 2))
        plt.ylim((0, 2))
        plt.grid()
        plt.plot(x_loc, u_all[i])

# %%
fig = plt.figure()

animation_1 = animation.FuncAnimation(fig, plotVelocity, interval=100)
plt.show()

# %%
