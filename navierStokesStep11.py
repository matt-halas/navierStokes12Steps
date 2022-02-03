# %%

import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt

nx = 51
x_loc = np.linspace(0, 2, nx)
ny = 51
y_loc = np.linspace(0, 2, ny)

dx = x_loc[1] - x_loc[0]
dy = y_loc[1] - y_loc[0]
dt = 0.001

nu = 0.01
rho = 1.23

def u_momentum(u, v, p, dt, dx, dy, nu, rho):
    u_int = u[1:-1, 1:-1]
    v_int = v[1:-1, 1:-1]
    un = np.zeros(np.shape(u))
    A = dt / dx**2 * (u[2:, 1:-1] - 2 * u_int + u[:-2, 1:-1]) \
        + dt / dy**2 * (u[1:-1, 2:] - 2 * u_int + u[1:-1, :-2])
    un[1:-1, 1:-1] = u_int \
        - u_int * dt / dx * (u_int - u[:-2, 1:-1]) \
        - v_int * dt / dy * (u_int - u[1:-1, :-2]) \
        - dt / (rho * 2 * dx) * (p[2:, 1:-1] - p[:-2, 1:-1]) \
        + nu * A
    un[:, -1] = 1
    un[:, 0] = 0
    un[(0, -1), :] = 0
    return un

def v_momentum(u, v, p, dt, dx, dy, nu, rho):
    u_int = u[1:-1, 1:-1]
    v_int = v[1:-1, 1:-1]
    vn = np.zeros(np.shape(u))
    A = dt / dx**2 * (v[2:, 1:-1] - 2 * v_int + v[:-2, 1:-1]) \
        + dt / dy**2 * (v[1:-1, 2:] - 2 * v_int + v[1:-1, :-2])
    vn[1:-1, 1:-1] = v_int \
        - u_int * dt / dx * (v_int - v[:-2, 1:-1]) \
        - v_int * dt / dy * (v_int - v[1:-1, :-2]) \
        - dt / (rho * 2 * dx) * (p[2:, 1:-1] - p[:-2, 1:-1]) \
        + nu * A
    vn[:, (0, -1)] = 0
    vn[(0, -1), :] = 0
    return vn

def pressure(u, v, p, dt, dx, dy, nu, rho, n_iterations):
    pn = np.zeros(np.shape(p))

    for i in range(n_iterations):
        p = pn.copy()

        A = 1 / dt * ((u[2:, 1:-1] - u[:-2, 1:-1] / (2 * dx)) \
            + (v[1:-1, 2:] - v[1:-1, :-2]) / (2 * dy)) \
            - ((u[2:, 1:-1] - u[:-2, 1:-1]) / (2 * dx))**2 \
            - 2 * ((u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dy) \
            * (v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dx)) \
            - ((v[1:-1, 2:] - v[1:-1, :-2]) / (2 * dy))**2
        pn[1:-1, 1:-1] = (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]) \
            + dx**2 * (p[1:-1, 2:] + p[1:-1, :-2])) \
            / 2 * (dx**2 + dy**2) \
            - rho * dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * A
        pn[(0, -1), :] = pn[(1, -2), :]
        pn[:, 0] = pn[:, 1]
        pn[:, -1] = 0
    
    return pn

# %%
p = np.zeros((nx, ny))
u = np.zeros((nx, ny))
v = np.zeros((nx, ny))

n_time_steps = 100
n_pres_iterations = 100

n_record = 50
record_interval = n_time_steps // n_record

p_all = np.zeros((n_record, nx, ny))
record_idx = 0

for i in range(n_time_steps):
    p = pressure(u, v, p, dt, dx, dy, nu, rho, n_pres_iterations)
    u = u_momentum(u, v, p, dt, dx, dy, nu, rho)
    v = v_momentum(u, v, p, dt, dx, dy, nu, rho)
    
    if record_interval != 0:
        if i % record_interval == 0:
            p_all[record_idx] = p
            record_idx += 1

# %%
def plotVelocity(i):
    if i >= len(p_all):
        plt.close()
    else:
        plt.clf()
        plt.xlim((0, 2))
        plt.ylim((0, 2))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.contourf(x_loc, y_loc, p_all[i].T, cmap='plasma',
            levels=25, vmin=-0.5, vmax=0.5)
        plt.colorbar()


if __name__ == '__main__':
    fig = plt.figure()
    animation_1 = animation.FuncAnimation(fig, plotVelocity, interval=50)
    animation_1.save('lidFlow.gif')