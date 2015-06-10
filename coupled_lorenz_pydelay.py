import numpy as np
import matplotlib.pyplot as plt
from pydelay import dde23


# define the equations
eqns = {
    'x1' : 'sigma*(y1 - x1)',
    'y1' : 'x1 * (rho1 - z1) - y1 + gamma11*pow(y1(t-delta11),2) + gamma12*pow(y2(t-delta12),2)',
    'z1' : 'x1 * y1 - beta*z1',
    'x2' : 'sigma*(y2 - x2)',
    'y2' : 'x2 * (rho2 - z2) - y2 + gamma21*pow(y1(t-delta21),2)',
    'z2' : 'x2 * y2 - beta*z2'
    }

#define the parameters
params = {
    'sigma': 10,
    'beta'  : 8/3.,
    'rho1' : 28,
    'rho2' : 28,
    'gamma11' : 0,
    'gamma12' : 0.1,
    'gamma21' : 0.1,
    'delta11' : 0,
    'delta12' : 45,
   	'delta21' : 75
    }

# Initialise the solver
dde = dde23(eqns=eqns, params=params)

# set the simulation parameters
# (solve from t=0 to t=1000 and limit the maximum step size to 1.0)
dde.set_sim_params(tfinal=1000, dtmax=0.01)

# set the history of to the constant function 0.5 (using a python lambda function)
histfunc = {
    'x1': lambda t: 0.,
    'y1' : lambda t: -0.05,
    'z1' : lambda t: 19,
    'x2' : lambda t: 0,
    'y2' : lambda t: 0.05,
    'z2' : lambda t: 19
    }

dde.hist_from_funcs(histfunc, 101)

# run the simulator
dde.run()

# Make a plot of x(t) vs x(t-tau):
# Sample the solution twice with a stepsize of dt=0.1:

# once in the interval [515, 1000]
dt = 0.001
t_from = 300*dt
t_to = (300 + 131072)*dt
sol1 = dde.sample(t_from, t_to, dt)
y1 = sol1['y1']
y2 = sol1['y2']
t = sol1['t']

dumps = np.zeros((t.shape[0], 3))
dumps[:,0] = t
dumps[:,1] = y1
dumps[:,2] = y2

np.savetxt("delay_lorenz.txt", dumps, fmt = "%.3f %.4f %.4f")

# f, ax = plt.subplots(2, sharex = True)
# ax[0].plot(t, y1, linewidth = 2, color = "#647C89")
# ax[1].plot(t, y2, linewidth = 2, color = "#647C89")
# ax[1].set_xlabel("time", fontsize = 25)
# ax[0].set_ylabel("y1", fontsize = 25)
# ax[1].set_ylabel("y2", fontsize = 25)
# plt.setp(ax[1].get_xticklabels(), fontsize=18)

# for i in range(2):
# 	plt.setp(ax[i].get_yticklabels(), fontsize=18)
# 	ax[i].set_xlim([0, t[-1]])
# plt.suptitle("Coupled delayed Lorenz system -- test case V in paper", fontsize = 38)

# plt.show()