import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import date, datetime
# from matplotlib.ticker import MultipleLocator, FuncFormatter
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
sys.path.append('/home/nikola/Work/multi-scale')
from src.data_class import DataField, load_enso_index

def plot_cycles(cycles, tit, fname = None):
    plt.figure(figsize=(15,10))
    plt.subplot(2, 1, 1)
    colors = ['k', '#7A1913'] + [np.random.rand(3,) for i in range(18,31,2)]
    widths = [1.2, 1.2] + [0.8 for i in range(18,31,2)]
    for cycle, col, wid in zip(cycles, colors, widths):
        plt.plot(cycle, color = col, linewidth = wid)
    plt.xlim([0, n34.data.shape[0]])
    plt.ylabel("annual, 5yr and QB", size = 20)
    plt.xticks(np.arange(0, n34.data.shape[0]+1, ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10))
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
    plt.xticks(np.arange(0, n34.data.shape[0]+1, ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10))
    plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)    
    plt.gca().spines['right'].set_visible(False)  
    plt.suptitle(tit, size = 25)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = 'tight')

def get_seasonal_indices(ts):
    # seasonal variance index
    m, var, d = ts.get_seasonality()
    var_month = np.power(var[:12], 2) # first twelve months time series must start in JANUARY
    v = (max(var_month) - min(var_month)) / max(var_month)
    ts.return_seasonality(m, var, d)

    # 2:1 phase synch
    import scipy.signal as ss
    analytic = ss.hilbert(ts.data)
    amp = np.sqrt(np.power(np.real(analytic), 2) + np.power(np.imag(analytic), 2))
    # phi_e = np.angle(analytic) % (2*np.pi)
    phi_e = 2*np.arctan2(np.imag(analytic), np.sqrt(np.power(np.real(analytic), 2) + np.power(np.imag(analytic), 2)) + np.real(analytic))
    phi_a = np.array([(2*np.pi/12) * t % (2*np.pi) for t in range(phi_e.shape[0])])
    xi = np.absolute(np.mean(np.exp(1j*(2*phi_e - phi_a))))

    # complex amplitude modulation
    psi = np.absolute(np.mean(amp*np.exp(1j*phi_a)))
    psi /= np.mean(amp)

    return v, xi, psi



FROM = 1950
TO = 2011
PERIODS = [12, 5*12] + [i for i in range(18,31,2)]

# ENSO REAL DATA -- Nino3.4 index
enso = load_enso_index("nino34raw.txt", '3.4', date(FROM, 1, 1), date(TO, 1, 1), anom = True)

## EMR
import scipy.io as sio
raw = sio.loadmat("Sergey-Nino34-ERM-linear-SST-20PC-L3-multiplicative-seasonal-std.mat")['N34s']
emr = DataField(data = raw[:, 0])
emr.create_time_array(date_from = date(1900, 1, 1), sampling = 'm')
emr.select_date(date(FROM, 1, 1), date(TO, 1, 1))

## CMIP5
model_names = ['CanESM2', 'CCSM4', 'CNRMCM5', 'CSIROmk360', 'GFDLCM3', 'GISSE2Hp1', 'GISSE2Hp2', 'GISSE2Hp3',
    'GISSE2Rp1', 'GISSE2Rp2', 'GISSE2Rp3', 'HadGem2ES', 'IPSL_CM5A_LR', 'MIROC5', 'MRICGCM3']
cmip = {}
for model in model_names:
    g = DataField()
    raw = np.loadtxt("N34_CMIP5/N34_%s.txt" % model)
    g.data = np.mean(raw, axis = 1)
    g.create_time_array(date(1861, 1, 1), sampling = 'm')
    g.select_date(date(FROM, 1, 1), date(TO, 1, 1))
    g.anomalise()
    cmip[model] = g

# CREATE N34 INDEX
n34 = load_enso_index("nino34raw.txt", '3.4', date(1900, 1, 1), date(2014, 1, 1), anom = False)
n34.get_seasonality(detrend = True, base_period = [date(1951,1,1), date(2000,12,1)])
n34.select_date(date(FROM, 1, 1), date(TO, 1, 1))
print get_seasonal_indices(enso)

# cycles = []
# for period in PERIODS:
#     enso.wavelet(period, 'm', save_wave = True)
#     cycles.append(enso.amplitude * np.cos(enso.phase))


# plot_cycles(cycles, tit = 'Nino3.4')#, fname = "plots/cycles/n34.png")
