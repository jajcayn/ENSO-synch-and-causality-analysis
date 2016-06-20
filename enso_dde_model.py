import numpy as np
from pydelay import dde23

eqn = {
    'h' : '-tanh(kappa * h(t-tau)) + b*cos(2*pi*t)'
}

params = {
    'kappa' : 100.,
    'tau' : 0.58,
    'b' : 1,
    'pi' : np.pi
}

dde = dde23(eqns = eqn, params = params)
dde.set_sim_params(tfinal=50000, dtmax=1.0)

histfunc = {
    'h' : lambda t: 1.        
}

dde.hist_from_funcs(histfunc, 100)
dde.run()

sol = dde.sample(20000, 20100, 0.01)
t = sol['t']
h = sol['h']

import matplotlib.pyplot as plt
plt.plot(t, h)
plt.show()