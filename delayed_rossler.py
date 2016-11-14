from pydelay import dde23
import numpy as np
# import matplotlib.pyplot as plt
import csv


delay_list = [0, 6, 30]
freqs_list = [[1.015, 0.985], [0.5, 2.515], [2.515, 0.5]]
names = ['1:1', '1:5', '5:1']
epsilon_list = [[0, 0.25]]


delay = 0
freqs = [2.398, 1.199] # annual, and biennal frequencies when considering monthly sampling
name = '1:2'
eps = [0, 0.25]

print("computing %s Rossler for epsilons in %.3f - %.3f" % (name, eps[0], eps[1]))

# coupled Rossler eqs
eqns = {
    'x1' : '-x3 - x2*omega1', ## unidirectional 
    # 'x1' : '-x3 - x2*omega1 + eps*(x4 - x1)', ## bidirectional 
    'x2' : 'omega1*x1 + a*x2',
    'x3' : 'b + x3*(x1 - c)',
    'x4' : '-x6 - x5*omega2 + eps*(x1(t - delta) - x4)',
    'x5' : 'omega2*x4 + a*x5',
    'x6' : 'b + x6*(x4 - c)'
    }

eps_list = np.linspace(eps[0], eps[1], 100)
dt = 0.157*2

rossler = {}
cnt = 0
for epsilon in eps_list:
    # parameters
    params = {
        'omega1': freqs[0],
        'omega2': freqs[1],
        'a': 0.15,
        'b': 0.20,
        'c': 10.0, 
        'eps' : epsilon,
        'delta' : delay * dt
        }


    # Initialise the solver
    dde = dde23(eqns=eqns, params=params)

    # set the simulation parameters
    # (solve from t=0 to t=1000 and limit the maximum step size to 1.0)
    dde.set_sim_params(tfinal=50000, dtmax=1.0)

    histfunc = {
        'x1': lambda t: 11.120979,
        'x2' : lambda t: 17.496796,
        'x3' : lambda t: 51.023544,
        'x4' : lambda t: 11.120979,
        'x5' : lambda t: 17.496796,
        'x6' : lambda t: 51.023544
        }

    dde.hist_from_funcs(histfunc, delay+5)

    # run the simulator
    dde.run()

    t_from = 1000*dt
    t_to = (131072 + 1000)*dt
    sol1 = dde.sample(t_from, t_to, dt)
    x1 = sol1['x1']
    x4 = sol1['x4']
    t = sol1['t']

    result = np.zeros((t.shape[0], 2))
    result[:, 0] = x1
    result[:, 1] = x4

    rossler[epsilon] = result
    cnt += 1
    print cnt, eps_list.shape[0]

fname = "conceptualRossler%smonthlysampling_100eps0-%.2f.dat" % (name, eps[1])
f = open(fname, 'w')
writer = csv.writer(f, lineterminator = "\n")
cnt = 1
for epsilon in eps_list:
    writer.writerow(["#Realization number %d" % cnt])
    writer.writerow(["#Delayed %s Rossler -- %d steps" % (name, delay)])
    writer.writerow(["#eps1 = %.4f" % epsilon])
    writer.writerow(["#eps2 = 0.0000"])
    writer.writerow(["#count = %d" % rossler[epsilon].shape[0]])
    for i in range(rossler[epsilon].shape[0]):
        row = ("%.6f   %.6f" % (rossler[epsilon][i, 0], rossler[epsilon][i,1]))
        writer.writerow([row])
    cnt += 1

f.close()




# f, ax = plt.subplots(2, sharex = True)
# ax[0].plot(t[120000:122000], x1[120000:122000], linewidth = 2, color = "#647C89")
# ax[1].plot(t[120000:122000], x4[120000:122000], linewidth = 2, color = "#647C89")
# ax[1].set_xlabel("time", fontsize = 25)
# ax[0].set_ylabel("x1", fontsize = 25)
# ax[1].set_ylabel("x4", fontsize = 25)
# plt.setp(ax[1].get_xticklabels(), fontsize=18)

# for i in range(2):
# 	plt.setp(ax[i].get_yticklabels(), fontsize=18)
# 	ax[i].set_xlim([t[120000], t[122000]])
# plt.suptitle("Coupled delayed Rossler -- epsilon %.3f // delay %d steps" % (EPS, DELAY), fontsize = 38)

# plt.show()
