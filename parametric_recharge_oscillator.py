import numpy as np
from scipy import integrate
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

year = 360. # 1 year as 360 days -- 12*30
DAMPED = False
SIGMA = 0.04
STATS = True

subtit = " -- modulation $\epsilon$ = $\pi / 2yr$"
# subtit = ""


def PRO_model(t, y, damped = False):

    if damped:
        lam = lam0 + eps * np.cos(omega_a * t)
        noise = np.random.normal(0, SIGMA)
    else:
        lam = eps * np.cos(omega_a * t)
        noise = 0.

    return [-lam * y[0] + omega_e * y[1] + noise, -omega_e * y[0]]


## params
twopi = 2*np.pi
# annual frequency
omega_a = twopi / year
# omega_a = 1. / (((2*np.pi) / 12) / 30)
# print omega_a
# enso period
omega_e = twopi / (3.75 * year)
# modulation
eps = np.pi / (2. * year)
# if damped, parameter of damping
lam0 = 1. / (0.4 * year)

mod = []

temp = []
y0 = [0., 1.35]

temp.append(y0[0])  
t0 = 0.

dt = 1.
if not DAMPED:
    time = 60
    ### SciPy integrator - R-K order 4
    r = integrate.ode(PRO_model).set_integrator('dopri5')
    r.set_initial_value(y0, t0)
    r.set_f_params(False)
    while r.successful and r.t < time * year:
        r.integrate(r.t + dt)
        temp.append(r.y[0])
        mod.append(eps * np.cos(omega_a * r.t))

    temp = np.array(temp[:-1])
    ts = []

    for i in range(time * 12):
        ts.append(np.mean(temp[i*30:(i+1)*30]))

else:
    time = 200
    ### Euler-Maruyama integrator
    t_tmp = t0
    y = np.zeros((int(time*year), 2))
    y[0, :] = y0
    for i in range(1, int(time*year)):
        y[i, 0] = y[i-1, 0] + dt * PRO_model(t_tmp, y[i-1, :], DAMPED)[0]
        y[i, 1] = y[i-1, 1] + dt * PRO_model(t_tmp, y[i-1, :], DAMPED)[1]
        t_tmp += dt

    temp = y[140*year:, 0]

    ts = []

    for i in range(60 * 12):
        ts.append(np.mean(temp[i*30:(i+1)*30]))


ts = np.array(ts)
# print np.var(ts)
# print np.std(ts)

hilb = np.imag(hilbert(ts))
if STATS:
    var = np.zeros((12,))
    amp = np.zeros_like(var)
    for i in range(12):
        var[i] = np.var(ts[i::12], ddof = 1)
        amp[i] = np.mean([np.sqrt(ts[j]**2 + hilb[j]**2) for j in range(i, ts.shape[0], 12)])

    

    print phase_diff
    plt.plot(phase_diff)
    plt.show()


plt.plot(ts, linewidth = 6., color = "#151618")
plt.plot(hilb, linewidth = 4, color = "#6B5F5C")
plt.axis([0, ts.shape[0], -4, 4])
plt.xticks(np.arange(0, ts.shape[0]+12, 60), np.arange(0, ts.shape[0]/12 + 1, 5))
plt.xlabel("time [years]", size = 25)
plt.ylabel("T [$^{\circ}$C]", size = 25)
if DAMPED:
    plt.title("PRO numerical integration (damped case) \n Euler-Maruyama method for SDE // noise = N(0, %.2f)" % SIGMA, size = 35)
else:
    plt.title("PRO numerical integration (neutral case) \n explicit Runge-Kutta of order 4 %s" % subtit, size = 35)

for it in (plt.gca().get_xticklabels() + plt.gca().get_yticklabels()):
    it.set_fontsize(22)

plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
# plt.close()
plt.show()

# plt.figure()
# plt.plot(mod)
# plt.xlim(0, len(mod))
# plt.show()

