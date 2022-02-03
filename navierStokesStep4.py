# %%
import numpy as np
import matplotlib.pyplot as plt

n_faces = 51
domain_length = 2
dx = 2 / (n_faces - 1)
x_loc = np.linspace(0, domain_length, n_faces)

dt = 0.001
n_time_steps = 4000
nu = 0.03
u = np.zeros(n_faces)
u[int(len(u)*0.25):int(len(u)*0.5)] = 2


def burgers_eqn(dt, dx, nu, u):
    un = np.zeros(len(u))
    un[1:-1] = u[1:-1] - u[1:-1] * dt / dx * (u[1:-1] - u[:-2]) \
        + nu * dt / dx**2 * (u[2:] - 2 * u[1:-1] + u[:-2])
    un[-1] = u[-1] - u[-1] * dt / dx * (u[-1] - u[-2]) \
        + nu * dt / dx**2 * (u[0] - 2 * u[-1] + u[-2])
    un[0] = un[-1]
    return un

n_record = 100
record_interval = n_time_steps // n_record

u_all = np.zeros((n_record, len(u)))

record_idx = 0
for i in range(n_time_steps):
    u = burgers_eqn(dt, dx, nu, u)
    if i % record_interval == 0:
        plt.plot(x_loc, u)
        u_all[record_idx] = u
        record_idx += 1

# %%
def plotVelocity(i):
    if i >= len(u_all):
        plt.close()
    else:
        plt.clf()
        plt.xlim((0, 2))
        plt.ylim((0, 2))
        plt.xlabel('x')
        plt.ylabel('velocity')
        plt.grid()
        plt.plot(x_loc, u_all[i])

# %%
import matplotlib.animation as animation

fig = plt.figure()

animation_1 = animation.FuncAnimation(fig, plotVelocity, interval=50)

animation_1.save('CFD/1DBurgersEqn.gif')
# %%
