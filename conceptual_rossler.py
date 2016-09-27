import numpy as np
import csv
import sys
sys.path.append("/Users/nikola/work-ui/multi-scale")
import src.mutual_information as MI
import src.wavelet_analysis as wvlt
from src.surrogates import get_single_FT_surrogate
import matplotlib.pyplot as plt


def read_rossler(fname):
    import csv
    r = {}
    f = open(fname, 'r')
    reader = csv.reader(f, lineterminator = "\n")
    eps = None
    for row in reader:
        if 'eps1' in row[0]:
            if eps is not None:
                r[eps] = np.array(r[eps])
            eps = float(row[0][8:])
            r[eps] = []
        if not "#" in row[0]:
            n = row[0].split()
            r[eps].append([float(n[0]), float(n[1])])
    r[eps] = np.array(r[eps])
    f.close()

    return r


fname = "conceptualRossler1:2monthlysampling_100eps0-0.25.dat"
# r = read_rossler(fname)
# print np.sort(r.keys())
# x, y = r[0.2374][:, 0], r[0.2374][:, 1]
# np.savetxt("temp.txt", r[0.2374][14276:18372, :], fmt = "%.8f")
a = np.loadtxt("temp.txt")
# # print a.shape
xts, yts = a[:2048, 0], a[:2048, 1]

# wave, _, _, _ = wvlt.continous_wavelet(xts, 1, True, wvlt.morlet, dj = 0, s0 = 24, j1 = 0, k0 = 6.)
# phase_x = np.arctan2(np.imag(wave), np.real(wave))[0, :]

# wave, _, _, _ = wvlt.continous_wavelet(yts, 1, True, wvlt.morlet, dj = 0, s0 = 12, j1 = 0, k0 = 6.)
# phase_y = np.arctan2(np.imag(wave), np.real(wave))[0, :]

xts -= np.mean(xts)
yts -= np.mean(yts)
xts /= np.std(xts)
yts /= np.std(yts)

import scipy.signal as ss

maxima = ss.argrelextrema(yts[:200], np.greater)[0][1::2]

print maxima

plt.plot(xts[:200], color = "#109EB2", linewidth = 2.3)
plt.plot(yts[:200] - 3.5, color = "#E67373", linewidth = 2.3)
# plt.plot(maxima, yts[maxima] - 3.5, 'o')
for i in range(maxima.shape[0]):
    plt.axvline(maxima[i], 0, 1, color = 'grey', linewidth = 1.2, linestyle = '-.')

# plt.xlim([-20, 220])
plt.savefig("phase-synch-ex.eps")




# cmi1knn = []
# for tau in range(1,7):
    
#     x, y, z = MI.get_time_series_condition([xts, yts], tau = tau, reversed = False, dim_of_condition = 3, eta = 3)
#     cmi1knn.append(MI.knn_cond_mutual_information(x, y, z, k = 64, dualtree = True))

# cmi = np.mean(cmi1knn)
# print cmi

# a = np.loadtxt("temp1.txt")

# t = np.greater(cmi, a)
# print np.sum(t)/500.

# import matplotlib.pyplot as plt
# plt.hist(a, bins = 50, fc = "grey", ec = 'grey')
# plt.axvline(cmi, 0, 1, color = "red", linewidth = 2)
# plt.savefig("surrogates.eps", bbox_inches = 'tight')

# surrs = []
# for num in range(500):
#     print num
#     x_t = get_single_FT_surrogate(xts)
#     cmi_temp = []
#     for tau in range(1,7):
#         x, y, z = MI.get_time_series_condition([x_t, yts], tau = tau, reversed = False, dim_of_condition = 3, eta = 3)
#         cmi_temp.append(MI.knn_cond_mutual_information(x, y, z, k = 64, dualtree = True))
#     surrs.append(np.mean(cmi_temp))
# # print surrs
# np.savetxt("temp1.txt", np.array(surrs))



# import matplotlib.pyplot as plt
# plt.plot(x[100:200], label = "x")
# plt.plot(y[100:200], label = "y")
# plt.legend()
# plt.show()

## x - biennal, y - annual ## frequencies as 1:2

# from multiprocessing import Pool


# def _pool_ros(args):
#     k, xts, yts = args
#     print k
#     cmi1bin = []
#     cmi2bin = []
#     cmi1knn = []
#     cmi2knn = []
#     for tau in range(1,31):
#         # first direction -- X -> Y | Y -- y is annual = 12 -->> eta = 12/4 = 3
#         x, y, z = MI.get_time_series_condition([xts, yts], tau = tau, reversed = False, dim_of_condition = 3, eta = 3)
#         cmi1bin.append(MI.cond_mutual_information(x, y, z, algorithm = 'EQQ2', bins = 8, log2 = False))
#         cmi1knn.append(MI.knn_cond_mutual_information(x, y, z, k = 32, dualtree = True))

#         # second direction - Y -> X | X -- x is biennal = 24 -->> eta = 24/4 = 6
#         x, y, z = MI.get_time_series_condition([xts, yts], tau = tau, reversed = True, dim_of_condition = 3, eta = 6)
#         cmi2bin.append(MI.cond_mutual_information(x, y, z, algorithm = 'EQQ2', bins = 8, log2 = False))
#         cmi2knn.append(MI.knn_cond_mutual_information(x, y, z, k = 32, dualtree = True))

#     CMI1bin = np.mean(cmi1bin)
#     CMI2bin = np.mean(cmi2bin)
#     CMI1knn = np.mean(cmi1knn)
#     CMI2knn = np.mean(cmi2knn)

#     return (k, CMI1bin, CMI2bin, CMI1knn, CMI2knn)

# pool = Pool(5)
# args = [(k, r[k][:, 0], r[k][:, 1]) for k in r.keys()]
# res = pool.map(_pool_ros, args)
# pool.close()
# pool.join()

# res = np.array(res)
# print res.shape
# np.savetxt(fname[:-4] + "-result.txt", res, fmt = "%.8f")