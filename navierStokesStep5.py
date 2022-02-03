# %%
import numpy as np
import matplotlib.animation as animation

nx = 41
x_loc = np.linspace(0, 2, nx)
ny = 41
y_loc = np.linspace(0, 2, ny)

dx = x_loc[1] - x_loc[0]
dy = y_loc[1] - y_loc[0]
dt = 0.001
c = 1

u = np.ones((nx, ny))

u[int(nx*0.25):int(nx*0.5), int(ny*0.25):int(ny*0.5)] = 2
# %%
import matplotlib.pyplot as plt

plt.contourf(x_loc, y_loc, u)

def step(dx, dy, dt, c, u):
    un = np.zeros_like(u)
    u_int = u[1:-1, 1:-1]
    un[1:-1, 1:-1] = u_int - c * dt / dx * (u_int - u[:-2, 1:-1]) \
        - c * dt / dy * (u_int - u[1:-1, :-2])
    un[:, (0, -1)] = 1
    un[(0, -1), :] = 1
    return un

n_time_steps = 500

n_record = 50
record_interval = n_time_steps // n_record

u_all = np.zeros((n_record, nx, ny))
record_idx = 0

for i in range(n_time_steps):
    u = step(dx, dy, dt, c, u)
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
        plt.contourf(x_loc, y_loc, u_all[i])

fig = plt.figure()

animation_1 = animation.FuncAnimation(fig, plotVelocity, interval=50)

animation_1.save('CFD/2DConv.gif')
# %%
