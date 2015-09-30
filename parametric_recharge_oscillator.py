import numpy as np
from scipy import integrate
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

year = 360. # 1 year as 360 days -- 12*30
DAMPED = True
SIGMA = 0.2
STATS = True
ENSO_PER = 4.75

subtit = " -- ENSO period $\omega_{e}$ = $%syr^{-1}$" % str(ENSO_PER)
# subtit = ""


def PRO_model(t, y, damped = False):

    if damped:
        lam = lam0 + eps * np.cos(omega_a * t)
        noise = np.random.normal(0, SIGMA) if SIGMA > 0. else 0.
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
omega_e = twopi / (ENSO_PER * year)
start_y = 40
end_y = 100
# modulation
eps = 2*np.pi / (2. * year)
# if damped, parameter of damping
lam0 = twopi / (0.4 * year)

mod = []

temp = []
y0 = [0., 1.35]

temp.append(y0[0])  
t0 = 0.

dt = 1.
if not DAMPED:
    time = 100
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

    for i in range(start_y*12, end_y * 12):
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
    print "variance -- ", np.var(ts, ddof = 1)
    print "std -- ", np.std(ts, ddof = 1)
    var = np.zeros((12,))
    amp = np.zeros_like(var)
    for i in range(12):
        var[i] = np.var(ts[i::12], ddof = 1)
        amp[i] = np.mean([np.sqrt(ts[j]**2 + hilb[j]**2) for j in range(i, ts.shape[0], 12)])

    phase_diff = np.zeros_like(ts)
    for t in range(ts.shape[0]):
        phi_e = np.arctan2(hilb[t], ts[t]) % (2*np.pi)
        phi_a = (2*np.pi/12) * t % (2*np.pi)
        phase_diff[t] = 2*phi_e - phi_a


    plt.figure(figsize=(20,16))
    if DAMPED:
        plt.suptitle("PRO numerical integration (damped case) \n Euler-Maruyama method for SDE // noise = N(0, %.2f)" % SIGMA, size = 35)
    else:
        plt.suptitle("PRO numerical integration (neutral case) \n explicit Runge-Kutta of order 4 %s" % subtit, size = 35)
    ax1 = plt.subplot2grid((2,2), (0,0), colspan = 2)
    ax1.plot(ts, linewidth = 6., color = "#151618")
    ax1.plot(hilb, linewidth = 4, color = "#6B5F5C")
    ax1.axis([0, ts.shape[0], -4, 4])
    ax1.set_xticks(np.arange(0, ts.shape[0]+12, 60))
    ax1.set_xticklabels(np.arange(start_y, start_y+ts.shape[0]/12 + 1, 5))
    ax1.set_xlabel("time [years]", size = 25)
    ax1.set_ylabel("T [$^{\circ}$C]", size = 25)
    for it in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        it.set_fontsize(22)
    ax1.xaxis.set_minor_locator(MultipleLocator(12))

    ax2 = plt.subplot2grid((2,2), (1, 0))
    ax2.plot(var, linewidth = 4, color = "#3A4864", label = "$\sigma^{2}_{m}$")
    ax2.plot(amp, linewidth = 4, color = "#95B05E", label = r"$\alpha_{m}$")
    ax2.axis([0, 11, 0, 3])
    ax2.set_xticks(np.arange(0, 12, 1))
    ax2.set_xticklabels(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
    ax2.set_ylabel(r"$\sigma^{2}_{m}$ [$^{\circ}C^{2}$] // $\alpha_{m}$ [$^{\circ}C$]", size = 25)
    ax2.legend()
    for it in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        it.set_fontsize(22)

    ax3 = plt.subplot2grid((2,2), (1, 1))
    ax3.hist(phase_diff, bins = 12, range = (0, 2*np.pi), normed = True, fc = "#6B5F5C", ec = "#151618")
    ax3.axis([0, 2*np.pi, 0, 0.35])
    ax3.set_xticks(np.arange(0, 3*np.pi, np.pi))
    ax3.set_xticklabels(["0", "$\pi$", "$2\pi$"])
    ax3.set_ylabel("probability", size = 25)
    ax3.set_xlabel("$\delta\phi_{1,2}$", size = 25)
    for it in (ax3.get_xticklabels() + ax3.get_yticklabels()):
        it.set_fontsize(22)

    if not DAMPED:
        # fname = "PROtest/PROneutral_modulation=pidiv2yr.png"
        fname = "PROtest/test.png"
    else:
        fname = "PROtest/test.png"

    plt.savefig(fname, bbox_inches = 'tight')
else:
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
    plt.show()
