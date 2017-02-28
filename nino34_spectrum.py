import numpy as np
import matplotlib.pyplot as plt
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
from matplotlib.ticker import MultipleLocator, FuncFormatter
import cPickle
import sys
import matplotlib.gridspec as gridspec
sys.path.append('/Users/nikola/work-ui/multi-scale')
# sys.path.append("/home/nikola/Work/multi-scale")
import src.wavelet_analysis as wvlt
import src.mutual_information as MI
from src.data_class import DataField, load_enso_index
from src.surrogates import SurrogateField
from scipy import signal


enso = load_enso_index("/Users/nikola/work-ui/data/nino34raw.txt", '3.4', date(1870, 1, 1), date(2016, 1, 1), anom = False)
enso_first = load_enso_index("/Users/nikola/work-ui/data/nino34raw.txt", '3.4', date(1870, 1, 1), date(1943, 1, 1), anom = False)
enso_second = load_enso_index("/Users/nikola/work-ui/data/nino34raw.txt", '3.4', date(1943, 1, 1), date(2016, 1, 1), anom = False)

WVLT_SPAN = [5, 96]
scales = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)

wvlt_power = np.zeros((scales.shape[0],))
wvlt_power_first = np.zeros_like(wvlt_power)
wvlt_power_second = np.zeros_like(wvlt_power)

for sc, i in zip(scales, range(wvlt_power.shape[0])):
    enso.wavelet(sc, period_unit = 'm', cut = 1, save_wave = True)
    enso_first.wavelet(sc, period_unit = 'm', cut = 1, save_wave = True)
    enso_second.wavelet(sc, period_unit = 'm', cut = 1, save_wave = True)

    wvlt_power[i] = np.sum(np.power(np.abs(enso.wave), 2)) / float(enso.wave.shape[0])
    wvlt_power_first[i] = np.sum(np.power(np.abs(enso_first.wave), 2)) / float(enso_first.wave.shape[0])
    wvlt_power_second[i] = np.sum(np.power(np.abs(enso_second.wave), 2)) / float(enso_second.wave.shape[0])


f, Pxx_spec = signal.welch(enso.data, 1./2.628e+6, 'blackman', 1024, scaling = 'spectrum')
f *= 2.628e+6#3.154e+7
f = 1./f



plt.figure()
# plt.plot(f, Pxx_spec)
plt.plot(scales, wvlt_power, color = '#242632', linewidth = 2., label = "FULL: 1870 - 2016")
plt.plot(scales, wvlt_power_first, linestyle = '--', color = '#F04F3B', linewidth = 1.8, label = "1st half: 1870 - 1943")
plt.plot(scales, wvlt_power_second, linestyle = ':', color = '#0FD86D', linewidth = 1.8, label = "2nd half: 1943 - 2016")
plt.xlim([scales[0], scales[-1]])
plt.xticks(scales[7::12], [0, 1, 2, 3, 4, 5, 6, 7, 8])
# print plt.xticks()
plt.xlabel("PERIOD [years]", size = 16)
plt.ylabel("WAVELET POWER", size = 16)
plt.legend(loc = "lower center")
# plt.gca().grid()
plt.gca().xaxis.set_major_locator(MultipleLocator(12))
# plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
plt.gca().xaxis.set_minor_locator(MultipleLocator(6))
for vl in range(6,scales[-1],6):
    plt.gca().axvline(vl, color = "#9B9B9B", linewidth = 0.5)
# plt.show()
plt.savefig("spectrum.eps")