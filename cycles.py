import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import date, datetime
# from matplotlib.ticker import MultipleLocator, FuncFormatter
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
sys.path.append('/home/nikola/Work/multi-scale')
from src.data_class import DataField, load_enso_index, load_station_data

def plot_cycles(cycles, tit, fname = None):
    plt.figure(figsize=(15,10))
    plt.subplot(2, 1, 1)
    np.random.seed(276)
    colors = ['#242632', '#F04F3B', '#873481'] + ['#1D43EF', '#1F7FF9', '#27ABE2', '#1FF2F9', '#1DEFBA', '#0FD86D', '#12F83B']#[np.random.rand(3,) for i in range(18,31,2)]
    widths = [2, 2, 2] + [1.4 for i in range(18,31,2)]
    for cycle, col, wid in zip(cycles, colors, widths):
        plt.plot(cycle, color = col, linewidth = wid)
    plt.xlim([0, enso.data.shape[0]])
    x = np.arange(enso.find_date_ndx(date(1982, 9, 1)), enso.find_date_ndx(date(1983, 5, 1)), 1)
    plt.fill_between(x, -6, 6, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
    x = np.arange(enso.find_date_ndx(date(1997, 7, 1)), enso.find_date_ndx(date(1998, 4, 1)), 1)
    plt.fill_between(x, -6, 6, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
    x = np.arange(enso.find_date_ndx(date(1972, 8, 1)), enso.find_date_ndx(date(1973, 3, 1)), 1)
    plt.fill_between(x, -6, 6, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
    x = np.arange(enso.find_date_ndx(date(1973, 9, 1)), enso.find_date_ndx(date(1974, 3, 1)), 1)
    plt.fill_between(x, -6, 6, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
    x = np.arange(enso.find_date_ndx(date(1975, 9, 1)), enso.find_date_ndx(date(1976, 3, 1)), 1)
    plt.fill_between(x, -6, 6, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
    x = np.arange(enso.find_date_ndx(date(1988, 8, 1)), enso.find_date_ndx(date(1989, 3, 1)), 1)
    plt.fill_between(x, -6, 6, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
    plt.ylabel("annual, 5yr and QB", size = 28)
    plt.xticks(np.arange(0, enso.data.shape[0]+1, ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10), rotation = 30)
    plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
    for tick in plt.gca().xaxis.get_major_ticks():
        tick.label.set_fontsize(19)
    for tick in plt.gca().yaxis.get_major_ticks():
        tick.label.set_fontsize(19) 
    plt.gca().spines['top'].set_visible(False)
    # plt.grid()
    plt.gca().spines['bottom'].set_visible(False)    
    plt.gca().spines['right'].set_visible(False) 
    plt.subplot(2, 1, 2)
    plt.plot(enso.data, color = '#242632', linewidth = 2.5)
    plt.axhline(-2, color = "#27ABE2", linewidth = 0.7)
    plt.axhline(2, color = "#27ABE2", linewidth = 0.7)
    plt.axhline(0, color = "#27ABE2", linewidth = 0.7)
    x = np.arange(enso.find_date_ndx(date(1982, 9, 1)), enso.find_date_ndx(date(1983, 5, 1)), 1)
    plt.fill_between(x, -3, 3, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.7)
    x = np.arange(enso.find_date_ndx(date(1997, 7, 1)), enso.find_date_ndx(date(1998, 4, 1)), 1)
    plt.fill_between(x, -3, 3, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.7)
    x = np.arange(enso.find_date_ndx(date(1972, 8, 1)), enso.find_date_ndx(date(1973, 3, 1)), 1)
    plt.fill_between(x, -3, 3, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.7)
    x = np.arange(enso.find_date_ndx(date(1973, 9, 1)), enso.find_date_ndx(date(1974, 3, 1)), 1)
    plt.fill_between(x, -3, 3, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.7)
    x = np.arange(enso.find_date_ndx(date(1975, 9, 1)), enso.find_date_ndx(date(1976, 3, 1)), 1)
    plt.fill_between(x, -3, 3, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.7)
    x = np.arange(enso.find_date_ndx(date(1988, 8, 1)), enso.find_date_ndx(date(1989, 3, 1)), 1)
    plt.fill_between(x, -3, 3, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.7)
    plt.xlim([0, enso.data.shape[0]])
    for tick in plt.gca().xaxis.get_major_ticks():
        tick.label.set_fontsize(19)
    for tick in plt.gca().yaxis.get_major_ticks():
        tick.label.set_fontsize(19) 
    plt.ylim([-3, 3])
    plt.gca().yaxis.set_minor_locator(MultipleLocator(0.5))
    plt.ylabel("Nino3.4 [SD]", size = 28)
    plt.xticks(np.arange(0, enso.data.shape[0]+1, ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10), rotation = 30)
    plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)    
    plt.gca().spines['right'].set_visible(False)  
    # plt.grid()
    plt.suptitle(tit, size = 25)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = 'tight')
    plt.close()

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
TO = 2000
PERIODS = [12, 5*12, 6*12] + [i for i in range(18,31,2)]
print PERIODS

# ENSO REAL DATA -- Nino3.4 index
# enso = load_enso_index("nino34raw.txt", '3.4', date(FROM-10, 1, 1), date(TO+10, 1, 1), anom = True)

# import scipy.io as sio
# a = sio.loadmat("Nino34-delay-quad-10PCsel-L3-seasonal-d:7mon-k:50.mat")['N34s']
# emr = DataField()
# emr.data = a[0, -1332:]
# emr.create_time_array(date(1900,1,1))
# emr.select_date(date(FROM, 1, 1), date(TO, 1, 1))
# emr.anomalise()
# emr.get_seasonality(detrend = True)

# plt.plot(emr.data)
# plt.show()

## EMR
# import scipy.io as sio
# raw = sio.loadmat("Sergey-Nino34-ERM-linear-SST-20PC-L3-multiplicative-seasonal-std.mat")['N34s']
# emr = DataField(data = raw[:, 0])
# emr.create_time_array(date_from = date(1900, 1, 1), sampling = 'm')
# emr.select_date(date(FROM, 1, 1), date(TO, 1, 1))

## CMIP5
# model_names = ['CanESM2', 'CCSM4', 'CNRMCM5', 'CSIROmk360', 'GFDLCM3', 'GISSE2Hp1', 'GISSE2Hp2', 'GISSE2Hp3',
#     'GISSE2Rp1', 'GISSE2Rp2', 'GISSE2Rp3', 'HadGem2ES', 'IPSL_CM5A_LR', 'MIROC5', 'MRICGCM3']
# cmip = {}
# # raw = sio.loadmat("S-Nino34-cond-noise-std.mat")['N34s']
# for model in model_names:
# # for i in range(20):
#     raw = np.loadtxt("N34_CMIP5/N34_%s.txt" % model)
#     for i in range(raw.shape[1]):
#     # g.data = np.mean(raw, axis = 1)
#         g = DataField()
#         g.data = raw[:, i]
#         g.create_time_array(date(1861, 1, 1), sampling = 'm')
#         # g.data = raw[:, i]
#         # g.create_time_array(date(1900, 1, 1), sampling = 'm')
#         g.select_date(date(FROM, 1, 1), date(TO, 1, 1))
#         g.anomalise()
#         cmip[model+str(i)] = g

# raw = np.loadtxt("N34_CMIP5/N34_CNRMCM5.txt")
# g = DataField()
# g.data = raw[:, 1]
# g.create_time_array(date(1861, 1, 1), sampling = 'm')
# g.select_date(date(FROM, 1, 1), date(TO, 1, 1))
# g.anomalise()


# CREATE N34 INDEX
# n34 = load_enso_index("nino34raw.txt", '3.4', date(1900, 1, 1), date(2014, 1, 1), anom = False)
# n34.get_seasonality(detrend = True, base_period = [date(1951,1,1), date(2000,12,1)])
# n34.select_date(date(FROM, 1, 1), date(TO, 1, 1))


# n34.data = enso.data.copy()
# vv, xxi, ppsi = get_seasonal_indices(enso)

# plt.figure(figsize=(15,10))
# for model in cmip:
# # for i in range(20):
#     # v, xi, psi = get_seasonal_indices(cmip[i])
#     # plt.subplot(1,2,1)
#     # plt.xlabel("ENSO seasonal variance", size = 20)
#     # plt.ylabel("2:1 phase synch", size = 20)
#     # plt.plot(v, xi, 'o', markersize = 8)
#     # plt.plot(vv, xxi, 's', markersize = 11, color = 'k')
#     # plt.subplot(1,2,2)
#     # plt.xlabel("ENSO seasonal variance", size = 20)
#     # plt.ylabel("Complex amplitude modulation", size = 20)
#     # plt.plot(v, psi, 'o', markersize = 8)
#     # plt.plot(vv, ppsi, 's', markersize = 11, color = 'k')
#     print model
#     cycles = []
#     for period in PERIODS:
#         cmip[model].wavelet(period, 'm', save_wave = True)
#         cycles.append(cmip[model].amplitude * np.cos(cmip[model].phase))
#     n34.data = cmip[model].data.copy()
#     plot_cycles(cycles, tit = "%s ts %d" % (model[:-1], int(model[-1])) , fname = "plots/cycles/%s-ts%d.png" % (model[:-1], int(model[-1])))

# # plt.legend()
# # plt.show()
# plt.savefig("plots/cycles/STATsynch.png")

# cycles = []
enso = load_station_data('../data/TG_STAID000027.txt', date(1775, 1, 1), date(2015, 1, 1), True)
# enso.get_monthly_data()
ndx = enso.select_date(date(1800, 1, 1), date(2005, 1, 1), apply_to_data = False)


enso.wavelet(8, 'y')
enso.phase = enso.phase[ndx]

# cont phase
for i in range(enso.phase.shape[0] - 1):
    if np.abs(enso.phase[i+1] - enso.phase[i]) > 1:
        enso.phase[i+1:] += 2*np.pi

diffs = np.diff(enso.phase)
diffs = diffs[diffs < 10.]
# plt.plot(diffs)
plt.hist(diffs, bins = 30, normed = True, fc = "#B0F8FF", ec = "#B0F8FF")
print diffs.min(), diffs.max()
plt.axvline(np.mean(diffs), 0, 1, color = '#17BEDD', linewidth = 2)
# plt.plot(diffs)
# plt.axvline(2*np.pi/(8*365.25), 0, 1, color = 'grey', linewidth = 1)
# plt.axvline(2*np.pi/(7*365.25), 0, 1, color = 'grey', linewidth = 1)
plt.xticks(2*np.pi/(np.arange(9,6,-0.5)*365.25), np.arange(9,5,-0.5))
plt.xlim([2*np.pi/(9.5*365.25), 2*np.pi/(6.5*365.25)])
plt.tick_params(axis='both', which='major', labelsize=15)
# plt.show()
plt.savefig("hist.eps", bbox_inches = 'tight')

# def get_equidistant_bins(num):
#     return np.array(np.linspace(-np.pi, np.pi, num+1))
# # for period in PERIODS:
# #     enso.wavelet(period, 'm', save_wave = True)
# #     cycles.append(enso.amplitude[ndx] * np.cos(enso.phase[ndx]))
# # enso.data = enso.data[ndx]
# # enso.time = enso.time[ndx]
# # plot_cycles(cycles, tit = '', fname = "nino34-cycles.eps")

# bins = get_equidistant_bins(8)

# bin_ndx = ((enso.phase[ndx] >= bins[0]) & (enso.phase[ndx] <= bins[1]))
# plt.figure(figsize=(20,8))
# d = enso.smoothing_running_avg(3)
# plt.plot(d[ndx], linewidth = 3.2, color = '#17BEDD')
# plt.plot(enso.phase[ndx] - 10, linewidth = 2.1, color = '#B0F8FF')
# plt.plot(enso.amplitude[ndx] - 12, linewidth = 2.1, color = '#FF9770')
# plt.plot(enso.amplitude[ndx] * np.cos(enso.phase[ndx]), linewidth = 4, color = '#FEF9F0')
# plt.fill_between(np.arange(d[ndx].shape[0]), enso.phase[ndx].min() - 10., (enso.amplitude[ndx] * np.cos(enso.phase[ndx])).max(), where = enso.phase[ndx] <= bins[1], facecolor = '#B0D2E0', edgecolor = '#B0D2E0', alpha = 0.4)
# plt.xlim([-10, enso.data[ndx].shape[0]+10])
# plt.savefig("wvlt-test.eps")
# plt.show()



# plot_cycles(cycles, tit = 'Nino3.4')#, fname = "plots/cycles/n34.png")

