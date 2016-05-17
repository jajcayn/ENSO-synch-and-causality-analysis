import numpy as np
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
from src.data_class import load_enso_index, DataField
from datetime import date
import matplotlib.pyplot as plt

n34 = load_enso_index('nino34raw.txt', 'n34', date(1870,1,1), date(2006,1,1), False)
n34.smoothing_running_avg(points = 24, cut_edges = False, use_to_data = True)
n34.get_seasonality()
plt.plot(n34.data, linewidth = 1.5, color = 'k', label = 'N3.4 index')
print 'NINO3.4'
print 'over 2: ', n34.data[n34.data > 2].shape
print 'under 2: ', n34.data[n34.data < -2.].shape

# models = ['CanESM2', 'CCSM4', 'CNRMCM5', 'CSIROmk360', 'GFDLCM3', 'GISSE2Hp1', 'GISSE2Hp2', 'GISSE2Hp3',
    # 'GISSE2Rp1', 'GISSE2Rp2', 'GISSE2Rp3', 'HadGem2ES', 'IPSL_CM5A_LR', 'MIROC5', 'MRICGCM3']
models = ['HadGem2ES', 'GISSE2Hp1', 'IPSL_CM5A_LR']
g = DataField()
corrs = []
for model in models:
    raw = np.loadtxt("CMIP5_NINO34/N34_%s.txt" % model)
    raw = raw[-1632:, :]
    g.time = n34.time.copy()
    g.data = np.mean(raw, axis = 1)
    g.smoothing_running_avg(points = 24, cut_edges = False, use_to_data = True)
    g.get_seasonality()

    print model
    print 'over 2: ', g.data[g.data > 2].shape
    print 'under 2:', g.data[g.data < -2].shape
    print 'corr:', np.corrcoef(g.data, n34.data)[0,1]
    corrs.append([np.corrcoef(g.data, n34.data)[0,1], model])

    plt.plot(g.data, linewidth = 0.8, label = "%s -- %.3f" % (model, np.corrcoef(g.data, n34.data)[0,1]))


print sorted(corrs, reverse = True)
_, _, y = n34.extract_day_month_year()
plt.xlim([0, n34.data.shape[0]])
t = np.arange(0, n34.data.shape[0], 60)
plt.axhline(2, 0, 1, color = 'r', linewidth = 2)
plt.axhline(-2, 0, 1, color = 'r', linewidth = 2)
plt.xticks(t, y[t], rotation = 45)
plt.legend()
plt.gca().grid()
plt.show()

