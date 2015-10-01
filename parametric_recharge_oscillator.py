import numpy as np
from scipy import integrate
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


class ENSO_PROmodel():

    def __init__(self, damped = False, daily = True, ENSOperiod = 3.75, modulation = 2, lambda0 = 0.4):

        self.damped = damped
        if daily:
            self.year = 360

        ## params
        twopi = 2*np.pi
        # annual frequency
        self.omega_a = twopi / self.year
        # enso period
        self.omega_e = twopi / (ENSOperiod * self.year)
        # modulation
        self.eps = 2*np.pi / (modulation * self.year)
        # if damped, parameter of damping
        self.lam0 = twopi / (lambda0 * self.year)


    def _PRO_model(self, t, y, damped = False, sigma = 0.04):

        if self.damped:
            lam = self.lam0 + self.eps * np.cos(self.omega_a * t)
            noise = np.random.normal(0, sigma) if sigma > 0. else 0.
        else:
            lam = self.eps * np.cos(self.omega_a * t)
            noise = 0.

        return [-lam * y[0] + self.omega_e * y[1] + noise, -self.omega_e * y[0]]


    def integrate_PROmodel(self, sigma = 0.04):
        
        # initial conditions
        temp = []
        y0 = [0., 1.35]
        temp.append(y0[0])  
        t0 = 0.
        dt = 1.

        if not self.damped:
            time = 60
            ### SciPy integrator - R-K order 4
            r = integrate.ode(self._PRO_model).set_integrator('dopri5')
            r.set_initial_value(y0, t0)
            r.set_f_params(False)
            while r.successful and r.t < time * self.year:
                r.integrate(r.t + dt)
                temp.append(r.y[0])

            temp = np.array(temp[:-1])
            ts = []

            for i in range(time * 12):
                ts.append(np.mean(temp[i*30:(i+1)*30]))

        else:
            time = 200
            ### Euler-Maruyama integrator
            t_tmp = t0
            y = np.zeros((int(time*self.year), 2))
            y[0, :] = y0
            for i in range(1, int(time*self.year)):
                y[i, 0] = y[i-1, 0] + dt * self._PRO_model(t_tmp, y[i-1, :], self.damped, sigma = sigma)[0]
                y[i, 1] = y[i-1, 1] + dt * self._PRO_model(t_tmp, y[i-1, :], self.damped, sigma = sigma)[1]
                t_tmp += dt

            temp = y[140*self.year:, 0]

            ts = []

            for i in range(60 * 12):
                ts.append(np.mean(temp[i*30:(i+1)*30]))

        self.time_series = np.array(ts)
        return np.array(ts)


year = 360. # 1 year as 360 days -- 12*30
DAMPED = False
STATS = True
subtit = ""
start_y = 0
end_y = 60

enso = ENSO_PROmodel()
enso.integrate_PROmodel()
ts = enso.time_series

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
