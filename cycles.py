"""
This scripts does various conditional analyses.
"""



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import date, datetime
# from matplotlib.ticker import MultipleLocator, FuncFormatter
from pyclits.geofield import DataField
from pyclits.data_loaders import load_enso_index
import scipy.stats as sts
import scipy.signal as ss


def get_frequency(phase, diff = True, window = 12):
    # if diff is True, simple differences are taken
    # else robust regression in windows

    # continuous phase
    ph = phase.copy()
    for t in range(ph.shape[0]-1):
        if np.abs(ph[t+1] - ph[t]) > 1:
            ph[t+1:] += 2*np.pi

    if diff:
        diffs = np.diff(ph)
        diffs = np.append(diffs, diffs[-1])
        diffs *= 12 / (2*np.pi)
    else:
        from sklearn import linear_model
        diffs = np.zeros_like(ph)
        for t in range(diffs.shape[0]):
            temp_ph = ph[max(t-window//2, 0) : min(t+window//2, diffs.shape[0])]
            model_r = linear_model.RANSACRegressor(linear_model.LinearRegression())
            x = np.arange(0, temp_ph.shape[0])
            x = x[:, np.newaxis]
            model_r.fit(x, temp_ph)
            coef = model_r.estimator_.coef_
            diffs[t] = coef * (12 / (2*np.pi)) 


    return diffs



def get_equidistant_bins(n):
    return np.linspace(-np.pi, np.pi, n+1)

def plot_cycles(cycles, phases = None, tit = '', fname = None, bins = None, bins_ts = None):
    plt.figure(figsize=(15,10))
    plt.subplot(2, 1, 1)
    np.random.seed(276)
    colors = ['#242632', '#F04F3B', '#873481'] + ['#1D43EF', '#1F7FF9', '#27ABE2', '#1FF2F9', '#1DEFBA', '#0FD86D', '#12F83B']#[np.random.rand(3,) for i in range(18,31,2)]
    widths = [2, 2, 2] + [1.4 for i in range(18,31,2)]
    if phases is None:
        for cycle, col, wid in zip(cycles, colors, widths):
            plt.plot(cycle, color = col, linewidth = wid)
    else:
        for cycle, col, wid, ph in zip(cycles, colors, widths, phases):
            plt.plot(cycle, color = col, linewidth = wid)
            plt.plot(ph, '--', color = col, linewidth = wid/1.5)

    if bins is not None:
        for line in bins:
            plt.axhline(line, xmin = 0, xmax = 1, color = 'k', linewidth = 0.7)
    if bins_ts is not None:
        to_plot = cycles[-1]
        for t in range(to_plot.shape[0]):
            if bins_ts[t] == 0:
                plt.plot(t, to_plot[t], 'x', color = 'C0')
            elif bins_ts[t] == 1:
                plt.plot(t, to_plot[t], 'o', color = 'C1')
            elif bins_ts[t] == 2:
                plt.plot(t, to_plot[t], 'v', color = 'C2')

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
    plt.ylabel("annual and 5 yr", size = 28)
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



FROM = 1910
TO = 2001
PERIODS = [12, 5*12] #+ [i for i in range(18,31,2)]
print PERIODS

# ENSO REAL DATA -- Nino3.4 index
enso = load_enso_index("nino34raw.txt", '3.4', date(FROM-10, 1, 1), date(TO+10, 1, 1), anom = False)
# enso.get_seasonality()
# enso = load_enso_index("nino34raw.txt", '3.4', date(1900, 1, 1), date(2011, 1, 1), anom = True)

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
# phases = []
# for period in PERIODS:
#     if period == PERIODS[-1]:
#         enso.wavelet(period, 'm', save_wave = True, cut = 36)#, cut_data = True, cut_time = True)
#     else:
#         enso.wavelet(period, 'm', save_wave = True, cut = 36)
#     phases.append(enso.phase)
#     cycles.append(enso.amplitude * np.cos(enso.phase))

enso.wavelet(60, period_unit = 'm', cut = 36)#, cut_data = True, cut_time = True)
lf = enso.amplitude * np.cos(enso.phase)
lf_ph = enso.phase.copy()
bins_ts = np.zeros_like(lf_ph)

# bins_tmp = np.linspace(lf_ph.min(), lf_ph.max(), 4)
# bins = []
# bins = list(np.linspace(bins_tmp[0], bins_tmp[1], 3))
# bins += list(np.linspace(bins_tmp[2], bins_tmp[3], 3))
# print bins
# for t in range(bins_ts.shape[0]):
#     if (lf_ph[t] < bins[1]) or (lf_ph[t] > bins[4]):
#         bins_ts[t] = 0
#     elif (lf_ph[t] > bins[1] and lf_ph[t] < bins[2]) or (lf_ph[t] > bins[3] and lf_ph[t] < bins[4]):
#         bins_ts[t] = 1
#     else:
#         bins_ts[t] = 2

bins = np.linspace(lf.min(), lf.max(), 4)
print bins
for t in range(bins_ts.shape[0]):
    if lf[t] < bins[1]:
        bins_ts[t] = 0
    elif (lf[t] > bins[1]) and (lf[t] < bins[2]):
        bins_ts[t] = 1
    else:
        bins_ts[t] = 2

# enso.wavelet(12, period_unit = 'm', cut = 36, cut_data = True, cut_time = True, regress_amp_to_data = False)
# ac = enso.amplitude * np.cos(enso.phase)
# # fit_x = np.vstack([ac, np.ones(ac.shape[0])]).T
# # m, c = np.linalg.lstsq(fit_x, enso.data)[0]
# # ac = m * ac + c

# plt.figure(figsize = (15,10))
# seasonal_cycle = np.zeros((12, 4))
# for b, lab, col, ls in zip(range(4), ['LF-', 'LF0', 'LF+', 'overall'], ['#1f77b4', '#673d00', '#c03634', 'k'], ['--', '--', '--', '-']):
#     for mon in range(1,13):
#         sel = enso.select_months(months = [mon], apply_to_data = False)
#         if b == 3:
#             seasonal_cycle[mon-1, b] = np.nanmean(enso.data[sel].copy())
#             plt.text(11.05, seasonal_cycle[-1, b], lab, horizontalalignment = 'left', verticalalignment = 'baseline', fontsize = 19)
#             # seasonal_cycle[mon-1, b] = np.nanmean(ac[sel].copy())
#         else:
#             data_temp = enso.data[sel].copy()
#             # data_temp = ac[sel].copy()
#             sel2 = (bins_ts[sel] == b)
#             seasonal_cycle[mon-1, b] = np.nanmean(data_temp[sel2])
#             plt.text(11.05, seasonal_cycle[-1, b], lab, horizontalalignment = 'left', verticalalignment = 'top', fontsize = 19, 
#                 fontdict = {'color' : col})

#     plt.plot(seasonal_cycle[:, b], ls, linewidth = 2.4, color = col)
#     # plt.text(11.05, seasonal_cycle[-1, b], lab, horizontalalignment = 'left', verticalalignment = 'center', fontsize = 14)
#     # plt.text(0.8, 0.75, "%.3f$\pm$%.3f" % (np.mean(ac[ndx]), np.std(ac[ndx], ddof = 1)), horizontalalignment = 'center', 
# #         verticalalignment = 'center', transform = plt.gca().transAxes, fontsize = 14)

# # ACminus = seasonal_cycle[:, 0].copy()
# # ACplus = seasonal_cycle[:, 2].copy()
# # ACzero = seasonal_cycle[:, 1].copy()

# # Cminus = []
# # Cplus = []
# # Czero = []
# time = []

# for year in range(enso.get_date_from_ndx(0).year, enso.get_date_from_ndx(-1).year+1):
#     tmp_ndx = enso.select_date(date(year, 1, 1), date(year+1, 1, 1), apply_to_data = False)
#     # print enso.data[tmp_ndx]
#     time.append(year)
#     if year == 1982:
#         plt.plot(enso.data[tmp_ndx], "o-", linewidth = 1.2, color = "#0C742C", markersize = 5, label = "NINO3.4: %d" % year)
#     elif year == 1983:
#         plt.plot(enso.data[tmp_ndx], "s-", linewidth = 1.2, color = "#0C742C", markersize = 5, label = "NINO3.4: %d" % year)
#     elif year == 1997:
#         plt.plot(enso.data[tmp_ndx], "o-", linewidth = 1.2, color = "#4B1A81", markersize = 5, label = "NINO3.4: %d" % year)
#     elif year == 1998:
#         plt.plot(enso.data[tmp_ndx], "s-", linewidth = 1.2, color = "#4B1A81", markersize = 5, label = "NINO3.4: %d" % year)
    # Cminus.append(sts.pearsonr(enso.data[tmp_ndx], ACminus)[0])
    # Cplus.append(sts.pearsonr(enso.data[tmp_ndx], ACplus)[0])
    # Czero.append(sts.pearsonr(enso.data[tmp_ndx], ACzero)[0])
    # Cminus.append(sts.pearsonr(ac[tmp_ndx], ACminus)[0])
    # Cplus.append(sts.pearsonr(ac[tmp_ndx], ACplus)[0])
    # Czero.append(sts.pearsonr(ac[tmp_ndx], ACzero)[0])

# Cminus = np.array(Cminus)
# Cplus = np.array(Cplus)
# Czero = np.array(Czero)
# time = np.array(time)

# fs = 1./3.154e+7

# f, Pxx_spec = ss.welch(Cminus, 1./3.154e+7, 'flattop', 1024, scaling='spectrum')
# f *= 3.154e+7
# fft = np.abs(np.fft.rfft(Cminus, axis = 0))
# freqs = np.fft.rfftfreq(Cminus.shape[0], d = 1./fs)
# freqs *= 3.154e+7
# plt.plot(freqs, 20*np.log10(fft), linewidth = 2., color = "#1f77b4", label = "spectrum C- full period")

# fft = np.abs(np.fft.rfft(Cminus[:Cminus.shape[0]//2], axis = 0))
# freqs = np.fft.rfftfreq(Cminus[:Cminus.shape[0]//2].shape[0], d = 1./fs)
# freqs *= 3.154e+7
# plt.plot(freqs, 20*np.log10(fft), ":", linewidth = 1.2, color = "#1f77b4", label = "spectrum C- 1st half")

# fft = np.abs(np.fft.rfft(Cminus[Cminus.shape[0]//2:], axis = 0))
# freqs = np.fft.rfftfreq(Cminus[Cminus.shape[0]//2:].shape[0], d = 1./fs)
# freqs *= 3.154e+7
# plt.plot(freqs, 20*np.log10(fft), "-.", linewidth = 1.2, color = "#1f77b4", label = "spectrum C- 2nd half")

# # f, Pxx_spec = ss.welch(Cplus, 1./3.154e+7, 'flattop', 1024, scaling='spectrum')
# # f *= 3.154e+7
# # plt.semilogy(f, np.sqrt(Pxx_spec), color = "#c03634", label = "spectrum C+")

# fft = np.abs(np.fft.rfft(Cplus, axis = 0))
# freqs = np.fft.rfftfreq(Cplus.shape[0], d = 1./fs)
# freqs *= 3.154e+7
# plt.plot(freqs, 20*np.log10(fft), linewidth = 2., color = "#c03634", label = "spectrum C+ full period")

# fft = np.abs(np.fft.rfft(Cplus[:Cplus.shape[0]//2], axis = 0))
# freqs = np.fft.rfftfreq(Cplus[:Cplus.shape[0]//2].shape[0], d = 1./fs)
# freqs *= 3.154e+7
# plt.plot(freqs, 20*np.log10(fft), ":", linewidth = 1.2, color = "#c03634", label = "spectrum C+ 1st half")

# fft = np.abs(np.fft.rfft(Cplus[Cplus.shape[0]//2:], axis = 0))
# freqs = np.fft.rfftfreq(Cplus[Cplus.shape[0]//2:].shape[0], d = 1./fs)
# freqs *= 3.154e+7
# plt.plot(freqs, 20*np.log10(fft), "-.", linewidth = 1.2, color = "#c03634", label = "spectrum C+ 2nd half")

# # f, Pxx_spec = ss.welch(Cplus - Cminus, 1./3.154e+7, 'flattop', 1024, scaling='spectrum')
# # f *= 3.154e+7
# # plt.semilogy(f, np.sqrt(Pxx_spec), ":",color = "#673d00", label = "spectrum C+ - C-")
# fft = np.abs(np.fft.rfft(Cplus - Cminus, axis = 0))
# freqs = np.fft.rfftfreq(Cminus.shape[0], d = 1./fs)
# freqs *= 3.154e+7
# plt.plot(freqs, 20*np.log10(fft), color = "#673d00", label = "spectrum C+ - C-")

# # ['#1f77b4', '#673d00', '#c03634']
# plt.plot(Cminus, linewidth = 2., label = "NINO34 | AC-", color = "#1f77b4")
# plt.plot(Cplus, linewidth = 2., label = "NINO34 | AC+", color = "#c03634")
# # plt.plot(Czero, linewidth = 2., label = "NINO3.4 | AC0", color = "#673d00")

# xlims = [-1, 1]

# x = np.array([np.where(time == 1982)[0][0], np.where(time == 1982)[0][0]+1])
# plt.fill_between(x, xlims[0], xlims[1], facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
# x = np.array([np.where(time == 1997)[0][0], np.where(time == 1997)[0][0]+1]) 
# plt.fill_between(x, xlims[0], xlims[1], facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
# x = np.array([np.where(time == 1972)[0][0], np.where(time == 1972)[0][0]+1]) 
# plt.fill_between(x, xlims[0], xlims[1], facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
# x = np.array([np.where(time == 1973)[0][0], np.where(time == 1973)[0][0]+1]) 
# plt.fill_between(x, xlims[0], xlims[1], facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
# x = np.array([np.where(time == 1975)[0][0], np.where(time == 1975)[0][0]+1]) 
# plt.fill_between(x, xlims[0], xlims[1], facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
# x = np.array([np.where(time == 1988)[0][0], np.where(time == 1988)[0][0]+1]) 
# plt.fill_between(x, xlims[0], xlims[1], facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)

# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.gca().spines['left'].set_visible(False)
# plt.xticks(np.arange(0,12), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], size = 15, rotation = 30)
# # plt.xticks(np.arange(0, Cminus.shape[0], 10), np.arange(enso.get_date_from_ndx(0).year, enso.get_date_from_ndx(-1).year, 10), size = 15, rotation = 30)
# plt.yticks(size = 15)
# # plt.ylabel("correlation NINO3.4 anom AC w/ AC anom cycles", size = 20)
# # # plt.xlabel("PERIOD [1/year]", size = 20)
# plt.ylabel("NINO3.4 [$^\circ$C]", size = 20)
# # # plt.ylabel("correlation NINO3.4 w/ AC cycles", size = 20)
# # # plt.ylim(xlims)
# plt.legend(prop = {'size' : 15}, frameon = False, loc = 8)
# # # plt.title("FFT spectrum, non-overlapping, RAW data", size = 25)
# plt.show()
# plt.savefig("plots/!resubmission/composites.eps", bbox_inches = 'tight')

print bins_ts.shape
# plot_cycles(cycles, phases, tit = '', fname = None, bins = bins, bins_ts = bins_ts)

# sst = load_ERSST_data('/Users/nikola/work-ui/data/ersstv4/', date(1900, 1, 1), date(2011,1,1), [-8, 8], [150, 270], anom = False)
# qb = [i for i in range(18,31,2)]
# qb_cyc = []
# for qb_per in qb:
#     sst.wavelet(qb_per, period_unit = 'm', cut = 36, regress_amp_to_data = True)
#     qb_cyc.append(sst.amplitude * np.cos(sst.phase))

# qb_cyc = np.array(qb_cyc)
# print qb_cyc.shape

# # sst.wavelet(24, period_unit = 'm', cut = 36, regress_amp_to_data = True)
# sst1y = np.mean(qb_cyc, axis = 0)
# # sst1y = sst.amplitude * np.cos(sst.phase)

# sst1ylf_minus = np.mean(sst1y[bins_ts == 2, ...], axis = 0)
# sst1ylf_0 = np.mean(sst1y[bins_ts == 1, ...], axis = 0)
# sst1ylf_plus = np.mean(sst1y[bins_ts == 0, ...], axis = 0)

# sst.quick_render(field_to_plot = sst1ylf_minus, symm = False, whole_world = False, tit = "LF-: SST QB cycle", 
#     fname = "plots/composites/SST-QB-LFminus-phase.png", cbar_label = "SST [C]")
# sst.quick_render(field_to_plot = sst1ylf_0, symm = False, whole_world = False, tit = "LF0: SST QB cycle", 
#     fname = "plots/composites/SST-QB-LF0-phase.png", cbar_label = "SST [C]")
# sst.quick_render(field_to_plot = sst1ylf_plus, symm = False, whole_world = False, tit = "LF+: SST QB cycle", 
#     fname = "plots/composites/SST-QB-LFplus-phase.png", cbar_label = "SST [C]")


enso.wavelet(12., period_unit = 'm', cut = 36, regress_amp_to_data = True)#, cut_data = True, cut_time = True)
# # print bins_ts.shape
# # enso.temporal_filter(cutoff = [18, 30], btype = 'bandpass', cut = 3, cut_time = True)
# # print enso.filtered_data.shape, enso.time.shape
# # from scipy.signal import hilbert
# # hilb = hilbert(enso.filtered_data)
# # enso.phase = np.arctan2(np.imag(hilb), np.real(hilb))

# ac = enso.amplitude * np.cos(enso.phase)
# ac_amp = enso.amplitude.copy()
# # ac0 = get_frequency(enso.phase, diff = True)
# # ac1 = get_frequency(enso.phase, diff = False, window = 6)
ac = get_frequency(enso.phase.copy(), diff = False, window = 12)
# ac2 = get_frequency(enso.phase.copy(), diff = False, window = 24)
# plt.plot(ac, "--", label = "12 window")
# plt.plot(ac2, label = "24 window")
# plt.plot(ac0, ":", label = "diffs")
# plt.plot(ac1, "-.",label = "6 window")
# plt.title("diffs: %.3f$\pm$%.3f | 6win: %.3f$\pm$%.3f | 12win: %.3f$\pm$%.3f | 24win: %.3f$\pm$%.3f" % (np.mean(ac0), 
#     np.std(ac0, ddof = 1), np.mean(ac1), np.std(ac1, ddof = 1), np.mean(ac), np.std(ac, ddof = 1), np.mean(ac2), 
#     np.std(ac2, ddof = 1)), size = 20)
# # plt.ylim([0, 2])
# # print np.mean(ac), np.mean(ac2)
# plt.xticks(np.arange(0, enso.data.shape[0], 10*12), np.arange(enso.get_date_from_ndx(0).year, enso.get_date_from_ndx(-1).year, 10))
# plt.axhline(1., 0, 1, color = 'k', linewidth = 1.6)
# plt.legend()
# plt.ylim([0.8, 1.2])
# x = np.arange(enso.find_date_ndx(date(1982, 9, 1)), enso.find_date_ndx(date(1983, 5, 1)), 1)
# plt.fill_between(x, 0.8, 1.2, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
# x = np.arange(enso.find_date_ndx(date(1997, 7, 1)), enso.find_date_ndx(date(1998, 4, 1)), 1)
# plt.fill_between(x, 0.8, 1.2, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
# x = np.arange(enso.find_date_ndx(date(1972, 8, 1)), enso.find_date_ndx(date(1973, 3, 1)), 1)
# plt.fill_between(x, 0.8, 1.2, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
# x = np.arange(enso.find_date_ndx(date(1973, 9, 1)), enso.find_date_ndx(date(1974, 3, 1)), 1)
# plt.fill_between(x, 0.8, 1.2, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
# x = np.arange(enso.find_date_ndx(date(1975, 9, 1)), enso.find_date_ndx(date(1976, 3, 1)), 1)
# plt.fill_between(x, 0.8, 1.2, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
# x = np.arange(enso.find_date_ndx(date(1988, 8, 1)), enso.find_date_ndx(date(1989, 3, 1)), 1)
# plt.fill_between(x, 0.8, 1.2, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
# plt.show()
# ac = np.nanmean(ac, axis = (1,2))
# ac_amp = get_frequency(enso.phase, diff = False, window = 24)
# ac_amp = np.nanmean(ac_amp, axis = (1,2))
# ac_amp = enso.amplitude.copy()
tits = ['LF-', 'LF0', 'LF+']
plt.figure(figsize = (14,8))
for b in range(3):
    plt.subplot(2,3,b+4)
    ndx = (bins_ts == b)
    weights = np.ones_like(ac[ndx])/float(len(ac[ndx]))
    plt.hist(ac[ndx], weights = weights, bins = 30, align = 'left', ec = '#666666', fc = '#666666')
    plt.text(0.8, 0.75, "%.3f$\pm$%.3f" % (np.mean(ac[ndx]), np.std(ac[ndx], ddof = 1)), horizontalalignment = 'center', 
        verticalalignment = 'center', transform = plt.gca().transAxes, fontsize = 14)
    # weights = np.ones_like(ac2[ndx])/float(len(ac2[ndx]))
    # plt.hist(ac2[ndx], weights = weights, bins = 30, align = 'left', alpha = 0.7, label = "24mon: %.3f$\pm$%.3f" % (np.mean(ac2[ndx]), np.std(ac2[ndx], ddof = 1)))
    # plt.legend(loc = 7, prop = {'size' : 9}, frameon = False)
    plt.xlabel("FREQUENCY [1/year]", size = 16)
    # plt.title(tits[b], size = 23)
    plt.gca().text(.8, 1., tits[b], horizontalalignment = 'right', verticalalignment = 'top', transform = plt.gca().transAxes, size = 24)
    plt.xlim([0.85, 1.15])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.axvline(1., 0, 1, color = 'k', linewidth = 1.2)
    plt.ylim([0, 0.1])
plt.subplot(2,3,2)
weights = np.ones_like(ac)/float(len(ac))
plt.hist(ac, weights = weights, bins = 30, align = 'left', ec = '#666666', fc = '#666666', label = "12mon: %.3f$\pm$%.3f" % (np.mean(ac), np.std(ac, ddof = 1)))
# weights = np.ones_like(ac2)/float(len(ac2))
# plt.hist(ac2, weights = weights, bins = 30, align = 'left', alpha = 0.7, label = "24mon: %.3f$\pm$%.3f" % (np.mean(ac2), np.std(ac2, ddof = 1)))
# plt.legend(loc = 7, prop = {'size' : 9}, frameon = False)
plt.text(0.8, 0.75, "%.3f$\pm$%.3f" % (np.mean(ac), np.std(ac, ddof = 1)), horizontalalignment = 'center', 
        verticalalignment = 'center', transform = plt.gca().transAxes, fontsize = 14)
# plt.suptitle("NINO3.4 AC frequencies", size = 25)
# plt.xlabel("PERIOD [1/year]", size = 16)
plt.xlim([0.85, 1.15])
plt.ylim([0, 0.1])
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.axvline(1., 0, 1, color = 'k', linewidth = 1.2)
# # plt.show()
plt.savefig("plots/!resubmission/ACfreq_LFcond.eps", bbox_inches = 'tight')

# print ac.shape, ac_amp.shape, lf_ph.shape

# enso.wavelet(24, period_unit = 'm', cut = 36, regress_amp_to_data = True)
# # # bc = get_frequency(enso.phase, diff = False, window = 12)
# bc = enso.amplitude * np.cos(enso.phase)
# # # bc = np.nanmean(bc, axis = (1,2))
# # # bc_amp = get_frequency(enso.phase, diff = False, window = 24)
# bc_amp = enso.amplitude.copy()
# # # bc_amp = np.nanmean(bc_amp, axis = (1,2))

# qbs = [i for i in range(18,31,2)]
# qb_cyc = []
# qb_amps = []
# for qb_per in qbs:
#     enso.wavelet(qb_per, period_unit = 'm', cut = 36, regress_amp_to_data = True)
# #     # qb_cyc.append(get_frequency(enso.phase, diff = False, window = 12))
#     qb_cyc.append(enso.amplitude * np.cos(enso.phase))
# #     # qb_amps.append(get_frequency(enso.phase, diff = False, window = 24))
#     qb_amps.append(enso.amplitude.copy())
# qb_cyc = np.array(qb_cyc)
# qb = np.mean(qb_cyc, axis = 0)
# # # qb = np.nanmean(qb, axis = (1,2))
# qb_amps = np.array(qb_amps)
# qb_amp = np.mean(qb_amps, axis = 0)
# # qb_amp = np.nanmean(qb_amp, axis = (1,2))


# bins = get_equidistant_bins(3)
# cond_means = np.zeros((3,6))

# for b in range(cond_means.shape[0]):
#     # ndx = ((lf_ph > bins[b]) & (lf_ph < bins[b+1]))
#     if b == 0:
#         ndx = (bins_ts == 0)
#     elif b== 1:
#         ndx = (bins_ts == 1)
#     else:
#         ndx = (bins_ts == 2)
#     cond_means[b, 0] = np.mean(ac[ndx])
#     cond_means[b, 1] = np.mean(ac_amp[ndx])
#     cond_means[b, 2] = np.mean(bc[ndx])
#     cond_means[b, 3] = np.mean(bc_amp[ndx])
#     cond_means[b, 4] = np.mean(qb[ndx])
#     cond_means[b, 5] = np.mean(qb_amp[ndx])

# diff = np.diff(bins)[0]
# # # diff = 1.
# plt.figure(figsize = (15,5))
# tits = ['AC', 'AC', 'BC', 'BC', 'QB', 'QB']
# for i in range(3):
#     # plt.subplot(1, 3, i+1)
#     # if i == 0:
#     #     # plt.ylabel("FREQ [1/year]", size = 18)
#     #     plt.ylabel("AMP [$^\circ$C]", size = 18)
#     # plt.bar(bins[:-1]+0.1*diff, cond_means[:, 2*i], width = 0.8*diff, align = 'edge', ec = 'k', fc = 'k')
#     # plt.title(tits[2*i], size = 20)
#     # plt.xticks(bins[:-1] + 0.5*diff, ['LF -', 'LF 0', 'LF +'], size = 17)
#     plt.subplot(1, 3, i+1)
#     if i == 0:
#         # plt.ylabel("FREQ [1/year]", size = 18)
#         plt.ylabel("AMPLITUDE [$^\circ$C]", size = 20)
#     plt.bar(bins[:-1]+0.1*diff, cond_means[:, 2*i+1], width = 0.8*diff, align = 'edge', ec = 'k', fc = 'k')
#     plt.title(tits[2*i+1], size = 25)
#     plt.gca().spines['top'].set_visible(False)
#     plt.gca().spines['right'].set_visible(False)
#     plt.gca().spines['left'].set_visible(False)
#     # plt.xlabel("LF phase", size = 18)
#     plt.ylim([0, 0.55])
#     plt.xticks(bins[:-1] + 0.5*diff, ['LF -', 'LF 0', 'LF +'], size = 20)
#     plt.yticks(size = 14)

    # plt.ylim([27, 28])

# plt.show()
# plt.savefig("N34anom_AMP_LFcond.eps", bbox_inches = 'tight')


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

