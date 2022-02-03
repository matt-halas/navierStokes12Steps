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
v = np.ones_like(u)

# %%
import matplotlib.pyplot as plt

plt.contourf(x_loc, y_loc, u)

def step(dx, dy, dt, c, u, v):
    un = np.zeros_like(u)
    vn = np.zeros_like(v)
    u_int = u[1:-1, 1:-1]
    v_int = v[1:-1, 1:-1]
    un[1:-1, 1:-1] = u_int - u_int * dt / dx * (u_int - u[:-2, 1:-1]) \
        - v_int * dt / dy * (u_int - u[1:-1, :-2])
    vn[1:-1, 1:-1] = v_int - u_int * dt / dx * (v_int - v[:-2, 1:-1]) \
        - v_int * dt / dy * (v_int - v[1:-1, :-2])
    un[:, (0, -1)] = 1
    un[(0, -1), :] = 1
    vn[:, (0, -1)] = 1
    vn[(0, -1), :] = 1
    return un, vn

n_time_steps = 500

n_record = 50
record_interval = n_time_steps // n_record

Vel_all = np.zeros((n_record, nx, ny))
u_all = np.zeros((n_record, nx, ny))
v_all = np.zeros((n_record, nx, ny))
record_idx = 0

quiver_step = 4
for i in range(n_time_steps):
    u, v = step(dx, dy, dt, c, u, v)
    if record_interval != 0:
        if i % record_interval == 0:
            Vel = np.sqrt(u**2 + v**2)
            plt.contourf(x_loc, y_loc, Vel)
            plt.quiver(x_loc[::quiver_step], y_loc[::quiver_step],
                u[::quiver_step, ::quiver_step], v[::quiver_step, ::quiver_step])
            Vel_all[record_idx] = Vel
            u_all[record_idx] = u
            v_all[record_idx] = v
            record_idx += 1

# %%
def plotVelocity(i):
    if i >= len(Vel_all):
        plt.close()
    else:
        plt.clf()
        plt.xlim((0, 2))
        plt.ylim((0, 2))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.contourf(x_loc, y_loc, Vel_all[i])
        plt.quiver(x_loc[::quiver_step], y_loc[::quiver_step],
            u_all[i, ::quiver_step, ::quiver_step],
            v_all[i, ::quiver_step, ::quiver_step])

if __name__ == '__main__':
    fig = plt.figure()
    animation_1 = animation.FuncAnimation(fig, plotVelocity, interval=200)
    plt.show()
    #animation_1.save('CFD/2DnonLinearConv.gif')
# %%
