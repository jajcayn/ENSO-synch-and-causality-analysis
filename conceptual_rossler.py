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
# fname = "ross-test.dat"
r = read_rossler(fname)
print np.sort(r.keys())
# x, y = r[0.2298][20000:52768, 0], r[0.2298][20000:52768, 1]
# np.savetxt("temp.txt", r[0.2449][20000:52768, :], fmt = "%.8f")
# a = np.loadtxt("temp.txt")
# # print a.shape
# xts, yts = a[:2048, 0], a[:2048, 1]


from multiprocessing import Pool


def _pool_ros(args):
    k, xts, yts = args
    print k
    cmi1bin = []
    cmi2bin = []
    cmi1knn = []
    cmi2knn = []
    for tau in range(1,31):
        # first direction -- X -> Y | Y -- y is biennal = 24 -->> eta = 24/4 = 3
        x, y, z = MI.get_time_series_condition([xts, yts], tau = tau, reversed = False, dim_of_condition = 3, eta = 6)
        cmi1bin.append(MI.cond_mutual_information(x, y, z, algorithm = 'EQQ2', bins = 8, log2 = False))
        cmi1knn.append(MI.knn_cond_mutual_information(x, y, z, k = 32, dualtree = True))

        # second direction - Y -> X | X -- x is annual = 24 -->> eta = 24/4 = 6
        x, y, z = MI.get_time_series_condition([xts, yts], tau = tau, reversed = True, dim_of_condition = 3, eta = 3)
        cmi2bin.append(MI.cond_mutual_information(x, y, z, algorithm = 'EQQ2', bins = 8, log2 = False))
        cmi2knn.append(MI.knn_cond_mutual_information(x, y, z, k = 32, dualtree = True))

    CMI1bin = np.mean(cmi1bin)
    CMI2bin = np.mean(cmi2bin)
    CMI1knn = np.mean(cmi1knn)
    CMI2knn = np.mean(cmi2knn)

    return (k, CMI1bin, CMI2bin, CMI1knn, CMI2knn)

pool = Pool(5)
args = [(k, r[k][20000:52768, 0], r[k][20000:52768, 1]) for k in r.keys()]
res = pool.map(_pool_ros, args)
pool.close()
pool.join()

res = np.array(res)
print res.shape
np.savetxt(fname[:-4] + "-result-32k.txt", res, fmt = "%.8f")