import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import date, datetime
# from matplotlib.ticker import MultipleLocator, FuncFormatter
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
from src.data_class import DataField, load_enso_index

def plot(cycles, tit, fname = None):
    plt.figure(figsize=(15,10))
    plt.subplot(2, 1, 1)
    colors = ['k', '#7A1913'] + [np.random.rand(3,) for i in range(18,31,2)]
    for cycle, col in zip(cycles, colors):
        plt.plot(cycle, color = col, linewidth = 0.8)
    plt.xlim([0, n34.data.shape[0]])
    plt.ylabel("annual, 5yr and QB", size = 20)
    plt.xticks(np.arange(0, n34.data.shape[0], ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10))
    plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)    
    plt.gca().spines['right'].set_visible(False) 
    plt.subplot(2, 1, 2)
    plt.plot(n34.data, 'k', linewidth = 2)
    plt.axhline(-2, color = "#01B8F1", linewidth = 1.2)
    plt.axhline(2, color = "#01B8F1", linewidth = 1.2)
    plt.axhline(0, color = "#01B8F1", linewidth = 1.2)
    plt.xlim([0, n34.data.shape[0]])
    plt.ylim([-3, 3])
    plt.gca().yaxis.set_minor_locator(MultipleLocator(0.5))
    plt.ylabel("Nino3.4 [SD]", size = 20)
    plt.xticks(np.arange(0, n34.data.shape[0], ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10))
    plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)    
    plt.gca().spines['right'].set_visible(False)  
    plt.suptitle(tit, size = 25)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = 'tight')



FROM = 1980
TO = 2000
PERIODS = [12, 5*12] + [i for i in range(18,31,2)]

enso = load_enso_index("nino34raw.txt", '3.4', date(FROM, 1, 1), date(TO, 1, 1), anom = False)
n34 = load_enso_index("nino34raw.txt", '3.4', date(1950, 1, 1), date(2001, 1, 1), anom = False)
n34.get_seasonality(base_period = [date(1951,1,1), date(2000,12,1)])
n34.select_date(date(FROM, 1, 1), date(TO, 1, 1))

cycles = []
for period in PERIODS:
    enso.wavelet(period, 'm', save_wave = True)
    cycles.append(enso.amplitude * np.cos(enso.phase))


plot(cycles, tit = 'NINO3.4 exp. data', fname = "plots/cycles/n34.png")

