import numpy as np
import matplotlib.pyplot as plt


def enso_dde_model(kappa, tau, b, nyears = 120, subsample_to_monthly = False):
    from pydelay import dde23

    eqn = {'h' : '-tanh(kappa * h(t-tau)) + b*cos(2*pi*t)'}
    params = {
        'kappa' : kappa,
        'tau' : tau,
        'b' : b,
        'pi' : np.pi
    }
    histfunc = {'h' : lambda t: 1.}
    
    dde = dde23(eqns = eqn, params = params)
    dde.set_sim_params(tfinal = 10000, dtmax = 0.01)
    dde.hist_from_funcs(histfunc, 100)
    dde.run()

    dt = 0.001

    sol = dde.sample(1000, 10000, dt) # 1 is year, int step is month
    t = sol['t']
    h = sol['h']

    t = t[:nyears/dt]
    h = h[:nyears/dt]

    if subsample_to_monthly:
        step = int((1/12.)/dt)
        t = t[::step]
        h = h[::step]

    return t, h


# t, h = enso_dde_model(50., 0.42, 1., nyears = 120, subsample_to_monthly = False)

# print h.shape

# plt.plot(t, h)
# plt.show()


