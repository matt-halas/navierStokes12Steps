import numpy as np
import matplotlib.pyplot as plt

n_faces = 1001
domain_length = 2
dx = 2 / (n_faces - 1)
x_loc = np.linspace(0, domain_length, n_faces)

dt = 0.001
n_time_steps = 200

u = np.ones(n_faces)
u[int(len(u)*0.25):int(len(u)*0.5)] = 2
du_wall = 0

def nonLinearStep(u, dt, dx, du_wall):
    u_diff = u[1:] - u[:-1]
    u_diff = np.append(du_wall, u_diff)
    u_new = u * (1 - dt / dx * u_diff)
    return u_new


for i in range(n_time_steps):
    u = nonLinearStep(u, dt, dx, du_wall)
    if i % 20 == 0:
        print(max(u * dt / dx))
        plt.plot(x_loc, u)