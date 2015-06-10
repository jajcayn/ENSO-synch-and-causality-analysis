import numpy as np
import matplotlib.pyplot as plt


def coupled_lorenz(t, y, v, n = 2):

    ui, vi, wi = [np.zeros((n,)) for i in range(3)]
    ui = y[0, :]
    vi = y[1, :]
    wi = y[2, :]

    for i in range(n):
        ui[i] = sigma * (vi[i] - ui[i])
        couple_sum = np.sum([gamma[iota,j] * np.power(v[j, t - delta[iota, j]], 2) for j in range(n) for iota in range(n)])
        vi[i] = ui[i] * (rho[i] - wi[i]) - vi[i] + couple_sum
        wi[i] = ui[i] * vi[i] - beta * wi[i]

    return [[ui[i], vi[i], wi[i]] for i in range(n)]




## params
N = 2

delta = np.zeros((N, N), dtype = np.int8)
delta[0, 0] = 0
delta[0, 1] = 45
delta[1, 0] = 75
delta[1, 1] = 0

gamma = np.zeros((N, N))
gamma[0, 0] = 0.
gamma[0, 1] = 0.1
gamma[1, 0] = 0.1
gamma[1, 1] = 0.

beta = 8/3.
sigma = 10.
rho = [28.0 for i in range(N)]

intSteps = 3000*50

t0 = 0.
dt = 0.0001
intStart = delta.max()
y0 = np.zeros((3, N))
y0[0, :] = 0.
y0[1, :] = [-0.05, 0.05]
y0[2, :] = 19. 

## integration
y = np.zeros((3, N, intSteps))

t_tmp = t0
for i in range(intStart):
    y[:, :, i] = y0 

## Euler-Maruyama
for i in range(intStart, intSteps):
    y[:, :, i] = y[:, :, i-1]
    result = coupled_lorenz(t_tmp, y[:, :, i-1], y[1, :, :], N)
    for dim in range(3):
        for iota in range(N):
            y[dim, iota, i] += dt * result[iota][dim]

    t_tmp += 1
    print t_tmp, y[:, 0, i]



## plot v - 2. dimension for both systems
f, ax = plt.subplots(N)
for i in range(N):
  ax[i].plot(y[1, i, :])

plt.show()
