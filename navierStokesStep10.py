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

p = np.zeros((nx, ny))

b = np.zeros((nx, ny))

b[int(1 / 4 * nx), int(1 / 4 * ny)] = 100
b[int(3 / 4 * nx), int(3 / 4 * ny)] = -100

# %%
def p_new(dx, dy, p, b):
    pn = np.zeros_like(p)
    pn[1:-1, 1:-1] = (dy**2 * (p[2:, 1:-1] + p[:-2, 1:-1]) \
        + dx**2 * (p[1:-1, 2:] + p[1:-1, :-2]) \
        - (b[1:-1, 1:-1] * dx**2 * dy**2)) \
        / (2 * (dx**2 + dy**2))
    pn[(0, -1), :] = 0
    pn[:, (0, -1)] = 0
    return pn

p_test = p_new(dx, dy, p, b)

n_time_steps = 100

n_record = 50
record_interval = n_time_steps // n_record

p_all = np.zeros((n_record, nx, ny))
record_idx = 0

for i in range(n_time_steps):
    p = p_new(dx, dy, p, b)
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
    animation_1.save('2DPoisson.gif')