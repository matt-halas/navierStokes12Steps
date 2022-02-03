# %%
import numpy as np
import matplotlib.animation as animation

nx = 81
x_loc = np.linspace(0, 2, nx)
ny = 81
y_loc = np.linspace(0, 2, ny)

dx = x_loc[1] - x_loc[0]
dy = y_loc[1] - y_loc[0]
nu = 0.05
sigma = 0.25
dt = sigma * dx * dy / nu

u = np.ones((nx, ny))

u[int(nx*0.25):int(nx*0.75), int(ny*0.25):int(ny*0.75)] = 2
# %%
import matplotlib.pyplot as plt

plt.contourf(x_loc, y_loc, u)

def step(dx, dy, dt, nu, u):
    un = np.zeros_like(u)
    u_int = u[1:-1, 1:-1]
    un[1:-1, 1:-1] = u_int \
        + nu * dt / dx**2 * (u[2:, 1:-1] - 2 * u_int + u[:-2, 1:-1]) \
        + nu * dt / dy**2 * (u[1:-1, 2:] - 2 * u_int + u[1:-1, :-2])
    un[:, (0, -1)] = 1
    un[(0, -1), :] = 1
    return un

n_time_steps = 500

n_record = 100
record_interval = n_time_steps // n_record

u_all = np.zeros((n_record, nx, ny))
record_idx = 0

for i in range(n_time_steps):
    u = step(dx, dy, dt, nu, u)
    if record_interval != 0:    
        if i % record_interval == 0:
            plt.contourf(x_loc, y_loc, u)
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
        plt.ylabel('y')
        plt.grid()
        plt.contourf(x_loc, y_loc, u_all[i], cmap='plasma',
            vmin=1, vmax=2, levels=100)
        plt.colorbar(ticks=np.linspace(1, 2, 6))

fig = plt.figure()

animation_1 = animation.FuncAnimation(fig, plotVelocity, interval=50)

animation_1.save('2DDiff.gif')
# %%
