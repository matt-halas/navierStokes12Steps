import numpy as np
import matplotlib.pyplot as plt

n_faces = 51
domain_length = 2
dx = 2 / (n_faces - 1)
x_loc = np.linspace(0, domain_length, n_faces)

dt = 0.0025
n_time_steps = 200

u = np.zeros(n_faces)
u[int(len(u)*0.25):int(len(u)*0.5)] = 2
u_wall = 1
nu = 0.3

def timeStep(u, nu, dt, dx):
    u_wall_left = np.interp(-dx, x_loc, u)
    u_wall_right = np.interp(2 + dx, x_loc, u)
    ui_plus = np.append(u, u_wall_right)
    ui_minus = np.append(u_wall_left, u)
    u_diff = ui_plus[1:] - 2 * u + ui_minus[:-1]
    return nu * dt / dx**2 * u_diff + u

for i in range(n_time_steps):
    u = timeStep(u, nu, dt, dx)
    if i % 10 == 0:
        plt.plot(x_loc, u)
