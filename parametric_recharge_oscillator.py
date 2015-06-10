import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

year = 360. # 1 year as 360 days -- 12*30
DAMPED = True
time = 60



def PRO_model(t, y, damped = False):

	if damped:
		lam = lam0 + eps * np.cos(omega_a * t)
		noise = np.random.normal(0, 0.04)
	else:
		lam = eps * np.cos(omega_a * t)
		noise = 0.

	return [-lam * y[0] + omega_e * y[1] + noise, -omega_e * y[0]]


## params
twopi = 2 * np.pi
# annual frequency
omega_a = twopi / year
# omega_a = 1. / (((2*np.pi) / 12) / 30)
# print omega_a
# enso period
omega_e = twopi / (3.75 * year)
# modulation
# eps = twopi/(2. * year)
eps = 2. / year
# if dampes, parameter of damping
lam0 = 0.4 / year

temp = []
y0 = [0., 1.35]

temp.append(y0[0])	
t0 = 0.

if not DAMPED:
    time += 40
    ## SciPy integrator
    print("starting integration")
    r = integrate.ode(PRO_model).set_integrator('dopri5')
    r.set_initial_value(y0, t0)
    r.set_f_params(False)
    while r.successful and r.t < (time * year):
    	r.integrate(r.t + 1)
    	temp.append(r.y[0])
    print("integration finished")
    temp = np.array(temp)
    temp = temp[40*year: ]
    ts = []
    for i in range((time-40) * 12):
        ts.append(np.mean(temp[i*30:(i+1)*30]))

else:
    ### Euler-Maruyama integrator
    time += 150
    dt = 1.
    t_tmp = t0
    y = np.zeros((int(time*year), 2))
    y[0, :] = y0
    for i in range(1, int(time*year)):
    	y[i, 0] = y[i-1, 0] + dt * PRO_model(t_tmp, y[i-1, :], True)[0]
    	y[i, 1] = y[i-1, 1] + dt * PRO_model(t_tmp, y[i-1, :], True)[1]
    	t_tmp += dt
    temp = y[150*year:, 0]
    ts = []
    for i in range((time-150) * 12):
    	ts.append(np.mean(temp[i*30:(i+1)*30]))


ts = np.array(ts)
fname = ("ensoPRO%s.txt" %('damped' if DAMPED else 'neutral'))

# np.savetxt(fname, ts, fmt = "%.3f")

plt.plot(ts, linewidth = 2., color = "#230503")
plt.axis([0, ts.shape[0], -4, 4])
plt.xticks(np.arange(0, ts.shape[0] + 5, 5*12), np.arange(0, ts.shape[0]/12 + 1, 5))
plt.xlabel("time [years]")
plt.ylabel("T [$^{\circ}$C]")
if DAMPED:
    plt.title("PRO numerical integration (damped case) -- Euler-Maruyama method for SDE // noise = N(0, 0.04)")
else:
    plt.title("PRO numerical integration (neutral case) -- Runge-Kutta of order 4")
plt.show()

