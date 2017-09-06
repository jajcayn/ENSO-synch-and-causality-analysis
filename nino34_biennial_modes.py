import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import date, datetime
# from matplotlib.ticker import MultipleLocator, FuncFormatter
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
sys.path.append('/home/nikola/Work/multi-scale')
from src.data_class import DataField, load_enso_index, load_station_data, load_ERSST_data
import scipy.stats as sts
import scipy.signal as ss

FROM = 1910
TO = 2001
PERIODS = [12, 5*12] #+ [i for i in range(18,31,2)]
print PERIODS

# ENSO REAL DATA -- Nino3.4 index
enso = load_enso_index("nino34raw.txt", '3.4', date(FROM-10, 1, 1), date(TO+10, 1, 1), anom = False)
mean, var, trend = enso.get_seasonality(detrend = False)
enso.return_seasonality(mean, var, trend)
# plt.plot(enso.data[36:-36], label = "raw SST")

# enso.wavelet(12., period_unit = 'm', cut = 36, regress_amp_to_data = True)#, cut_data = True, cut_time = True)
enso.temporal_filter(cutoff = [18, 30], btype = 'bandpass', cut = 3)#, cut_time = True, cut_data = True)
lf = enso.filtered_data.copy()
# hilb = ss.hilbert(enso.filtered_data)
# enso.phase = np.arctan2(np.imag(hilb), np.real(hilb))
# enso.amplitude = np.sqrt(np.power(np.real(hilb),2) + np.power(np.imag(hilb),2))
# oscill = enso.amplitude * np.cos(enso.phase)

bins_ts = np.zeros_like(lf)
bins = np.linspace(lf.min(), lf.max(), 4)
print bins
for t in range(bins_ts.shape[0]):
    if lf[t] < bins[1]:
        bins_ts[t] = 0
    elif (lf[t] > bins[1]) and (lf[t] < bins[2]):
        bins_ts[t] = 1
    else:
        bins_ts[t] = 2


enso.wavelet(12*5, period_unit = 'm', cut = 36, regress_amp_to_data = True, cut_data = True, cut_time = True)
# hilb = ss.hilbert(enso.filtered_data)
# enso.phase = np.arctan2(np.imag(hilb), np.real(hilb))
# enso.amplitude = np.sqrt(np.power(np.real(hilb),2) + np.power(np.imag(hilb),2))
lf = enso.amplitude * np.cos(enso.phase)

bins_ts = np.zeros_like(lf)
bins = np.linspace(lf.min(), lf.max(), 4)
print bins
for t in range(bins_ts.shape[0]):
    if lf[t] < bins[1]:
        bins_ts[t] = 0
    elif (lf[t] > bins[1]) and (lf[t] < bins[2]):
        bins_ts[t] = 1
    else:
        bins_ts[t] = 2

# plt.plot(bins_ts)
# plt.show()

seasonal_cycle = np.zeros((12, 7))
for b, lab, col, ls in zip(range(4), ['QB-', 'QB0', 'QB+', 'overall'], ['#1f77b4', '#673d00', '#c03634', 'k'], ['--', '--', '--', '-']):
    for mon in range(1,13):
        sel = enso.select_months(months = [mon], apply_to_data = False)
        if b == 3:
            seasonal_cycle[mon-1, b] = np.nanmean(enso.data[sel].copy())
            plt.text(11.05, seasonal_cycle[-1, b], lab, horizontalalignment = 'left', verticalalignment = 'baseline', fontsize = 19)
            # seasonal_cycle[mon-1, b] = np.nanmean(ac[sel].copy())
        else:
            data_temp = enso.data[sel].copy()
            # data_temp = ac[sel].copy()
            sel2 = (bins_ts[sel] == b)
            seasonal_cycle[mon-1, b] = np.nanmean(data_temp[sel2])
            plt.text(11.05, seasonal_cycle[-1, b], lab, horizontalalignment = 'left', verticalalignment = 'top', fontsize = 19, 
                fontdict = {'color' : col})

    plt.plot(seasonal_cycle[:, b], ls, linewidth = 2.4, color = col)
    # plt.text(11.05, seasonal_cycle[-1, b], lab, horizontalalignment = 'left', verticalalignment = 'center', fontsize = 14)
    # plt.text(0.8, 0.75, "%.3f$\pm$%.3f" % (np.mean(ac[ndx]), np.std(ac[ndx], ddof = 1)), horizontalalignment = 'center', 
        # verticalalignment = 'center', transform = plt.gca().transAxes, fontsize = 14)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xticks(np.arange(0,12), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], size = 15, rotation = 30)
plt.ylabel("NINO3.4 [$^\circ$C]", size = 20)
plt.show()


