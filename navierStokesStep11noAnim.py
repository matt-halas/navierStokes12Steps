# %%
import numpy as np
import matplotlib.pyplot as plt

def momentum(u, v, p, dt, dx, dy, nu, rho):
    u_int = u[1:-1, 1:-1].copy()
    v_int = v[1:-1, 1:-1].copy()
    un = np.zeros(np.shape(u))
    A = dt / dx**2 * (u[1:-1, 2:] - 2 * u_int + u[1:-1, :-2]) \
        + dt / dy**2 * (u[2:, 1:-1] - 2 * u_int + u[:-2, 1:-1])
    un[1:-1, 1:-1] = u_int \
        - u_int * dt / dx * (u_int - u[1:-1, :-2]) \
        - v_int * dt / dy * (u_int - u[:-2, 1:-1]) \
        - dt / (rho * 2 * dx) * (p[1:-1, 2:] - p[1:-1, :-2]) \
        + nu * A
    un[-1, :] = 1
    un[0, :] = 0
    un[:, (0, -1)] = 0

    vn = np.zeros(np.shape(u))
    A = dt / dx**2 * (v[1:-1, 2:] - 2 * v_int + v[1:-1, :-2]) \
        + dt / dy**2 * (v[2:, 1:-1] - 2 * v_int + v[:-2, 1:-1])
    vn[1:-1, 1:-1] = v_int \
        - u_int * dt / dx * (v_int - v[1:-1, :-2]) \
        - v_int * dt / dy * (v_int - v[:-2, 1:-1]) \
        - dt / (rho * 2 * dy) * (p[2:, 1:-1] - p[:-2, 1:-1]) \
        + nu * A
    vn[:, (0, -1)] = 0
    vn[(0, -1), :] = 0

    return un, vn

def pressure(u, v, p, dt, dx, dy, nu, rho, n_iterations):
    for i in range(n_iterations):
        pn = p.copy()
        du_dx = (u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dx)
        dv_dy = (v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dy)
        du_dy = (u[2:, 1:-1] - u[:-2, 1:-1]) / (2 * dy)
        dv_dx = (v[1:-1, 2:] - v[1:-1, :-2]) / (2 * dx)
        A = (du_dx + dv_dy) / dt - du_dx**2 - 2 * du_dx * du_dy \
            - dv_dy**2
        p[1:-1, 1:-1] = (dy**2 * (pn[1:-1, 2:] + pn[1:-1, :-2]) \
            + dx**2 * (pn[2:, 1:-1] + pn[:-2, 1:-1])) \
            / (2 * (dx**2 + dy**2)) \
            - rho * dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * A
        p[-1, :] = 0 # pressure at y=2 is 0
        p[0, :] = p[1, :] # dp/dy is zero at y=0
        p[:, 0] = p[:, 1] # dp/dx is zero at x=0
        p[:, -1] = p[:, -2] # dp/dx is zero at x=2
    return p

# %%
nx = 81
x_max = 2
x_loc = np.linspace(0, x_max, nx)
ny = 81
y_max = 2
y_loc = np.linspace(0, y_max, ny)
dx = x_max / (nx - 1)
dy = y_max / (ny - 1)

X, Y = np.meshgrid(x_loc, y_loc)

dt = 0.001

nu = 0.1
rho = 0.1

p = np.zeros((ny, nx))
u = np.zeros((ny, nx))
v = np.zeros((ny, nx))

n_time_steps = 20
n_pres_iterations = 500

for i in range(n_time_steps):
    p = pressure(u, v, p, dt, dx, dy, nu, rho, n_pres_iterations)
    u, v = momentum(u, v, p, dt, dx, dy, nu, rho)

# %%

plt.figure()
plt.contourf(x_loc, y_loc, p, cmap='plasma')
plt.title(f"Pressure dist at {n_time_steps * dt} seconds")
plt.colorbar()
plt.streamplot(x_loc, y_loc, u, v, density=1)
plt.xlim(0, 2)
plt.ylim(0, 2)
# %%
