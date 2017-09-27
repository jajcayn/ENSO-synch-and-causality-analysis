"""
This scripts implements various analyses of biennial modes of ENSO
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import date, datetime
# from matplotlib.ticker import MultipleLocator, FuncFormatter
from pyclits.data_class import DataField
from pyclits.data_loaders import load_enso_index
import scipy.stats as sts
import scipy.signal as ss

FROM = 1910
TO = 2001
PERIODS = [12, 5*12] #+ [i for i in range(18,31,2)]
print PERIODS

# ENSO REAL DATA -- Nino3.4 index
enso = load_enso_index("nino34raw.txt", '3.4', date(FROM-10, 1, 1), date(TO+10, 1, 1), anom = False)
# mean, var, trend = enso.get_seasonality(detrend = False)
# enso.return_seasonality(mean, var, trend)
# plt.plot(enso.data[36:-36], label = "raw SST")

# enso.wavelet(12., period_unit = 'm', cut = 36, regress_amp_to_data = True)#, cut_data = True, cut_time = True)
enso.temporal_filter(cutoff = [18, 30], btype = 'bandpass', cut = 3)#, cut_time = True, cut_data = True)
qb = enso.filtered_data.copy()
# enso.anomalise()

bins_ts2 = np.zeros_like(qb)
bins = np.linspace(qb.min(), qb.max(), 4)
for t in range(bins_ts2.shape[0]):
    if qb[t] < bins[1]:
        bins_ts2[t] = 0
    elif (qb[t] > bins[1]) and (qb[t] < bins[2]):
        bins_ts2[t] = 1
    else:
        bins_ts2[t] = 2


enso.wavelet(60, period_unit = 'm', cut = 36, cut_data = True, cut_time = True)
lf = enso.amplitude * np.cos(enso.phase)
lf_ph = enso.phase.copy()
bins_ts = np.zeros_like(lf_ph)

bins = np.linspace(lf.min(), lf.max(), 4)
print bins
for t in range(bins_ts.shape[0]):
    if lf[t] < bins[1]:
        bins_ts[t] = 0
    elif (lf[t] > bins[1]) and (lf[t] < bins[2]):
        bins_ts[t] = 1
    else:
        bins_ts[t] = 2

qb -= np.mean(qb)

qbplus = qb * (qb > 0)
qbminus = qb * (qb < 0)


# plt.plot(enso.data, label = "LF")
# for t in range(bins_ts.shape[0]):
#     if bins_ts[t] == 0:
#         plt.plot(t, enso.data[t], 'x', color = 'blue', markersize = 4)
#     elif bins_ts[t] == 1:
#         plt.plot(t, enso.data[t], 'x', color = 'gray', markersize = 4)
#     elif bins_ts[t] == 2:
#         plt.plot(t, enso.data[t], 'x', color = 'red', markersize = 4)

# plt.plot(enso.data - 8, label = "QB")
# for t in range(bins_ts2.shape[0]):
#     if bins_ts2[t] == 0:
#         plt.plot(t, enso.data[t] - 8., 'o', color = 'blue', markersize = 4)
#     elif bins_ts2[t] == 1:
#         plt.plot(t, enso.data[t] - 8., 'o', color = 'gray', markersize = 4)
#     elif bins_ts2[t] == 2:
#         plt.plot(t, enso.data[t] - 8., 'o', color = 'red', markersize = 4)

# plt.xticks(np.arange(0, enso.time.shape[0], 120), np.arange(enso.get_date_from_ndx(0).year, 
#     enso.get_date_from_ndx(-1).year, 10), rotation = 30)
# plt.legend()
# plt.ylim([15,30])
# plt.yticks([], [])


# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.gca().spines['left'].set_visible(False)
# plt.show()

# plt.plot(qbplus, label = "QB+", color = "#c03634", linewidth = 2)
# plt.plot(qbminus, label = "QB-", color = "#1f77b4", linewidth = 2)
# plt.ylabel("QB [$^\circ$C]", size = 20)
# plt.yticks(size = 15)
# plt.xticks(np.arange(0, enso.time.shape[0], 10*12), np.arange(enso.get_date_from_ndx(0).year, 
#     enso.get_date_from_ndx(-1).year, 10), rotation = 30, size = 15)
# plt.legend()
# plt.show()
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.gca().spines['left'].set_visible(False)
# plt.xticks(np.arange(0,12), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], size = 15, rotation = 30)
# plt.ylabel("QB [$^\circ$C]", size = 20)
# plt.yticks(size = 15)
# plt.legend()
# plt.show()

# plt.plot(enso.data, label = "NINO3.4 anom", linewidth = 1.1, color = 'k')
# plt.plot(qb, label = "QB [18-30bandpass] NINO3.4", linewidth = 2., color = 'C1')
# plt.xticks(np.arange(0, enso.time.shape[0], 10*12), np.arange(enso.get_date_from_ndx(0).year, 
#     enso.get_date_from_ndx(-1).year, 10), rotation = 30)
# plt.legend()
# plt.show()
# hilb = ss.hilbert(enso.filtered_data)
# enso.phase = np.arctan2(np.imag(hilb), np.real(hilb))
# enso.amplitude = np.sqrt(np.power(np.real(hilb),2) + np.power(np.imag(hilb),2))
# oscill = enso.amplitude * np.cos(enso.phase)

# plt.plot(bins_ts)
# plt.plot(bins_ts2)
# plt.show()
# plt.subplot(211)
seasonal_cycle = np.zeros((12, 4))
# # labs = ['LF-', 'LF0', 'LF+', 'overall']
# labs = ['LF-', 'LF0', 'LF+', 'QB-', 'QB0', 'QB+', 'overall']
# # cols = ['#1f77b4', '#673d00', '#c03634', 'k']
cols = ['#1f77b4', '#673d00', '#c03634', '#1f77b4', '#673d00', '#c03634', 'k']
# # lss = ['--', '--', '--', '-']
# lss = ['--', '--', '--', ':', ':', ':', '-']
plt.figure()
# for b, lab, col, ls in zip(range(4), labs, cols, lss):
qbp_lfp = np.logical_and(bins_ts2 == 2, bins_ts == 2)
qbp_lfm = np.logical_and(bins_ts2 == 2, bins_ts == 0)
qbm_lfp = np.logical_and(bins_ts2 == 0, bins_ts == 2)
qbm_lfm = np.logical_and(bins_ts2 == 0, bins_ts == 0)
qbp_lfp_list = [] 
qbp_lfm_list = [] 
qbm_lfp_list = [] 
qbm_lfm_list = [] 
for mon in range(1,13):
    sel = enso.select_months(months = [mon], apply_to_data = False)
    data_temp = enso.data[sel].copy()
    qbp_lfp_list.append(np.sum(qbp_lfp[sel]))
    qbp_lfm_list.append(np.sum(qbp_lfm[sel]))
    qbm_lfp_list.append(np.sum(qbm_lfp[sel]))
    qbm_lfm_list.append(np.sum(qbm_lfm[sel]))
    seasonal_cycle[mon-1, 0] = np.nanmean(data_temp[qbp_lfp[sel]])
    seasonal_cycle[mon-1, 1] = np.nanmean(data_temp[qbp_lfm[sel]])
    seasonal_cycle[mon-1, 2] = np.nanmean(data_temp[qbm_lfp[sel]])
    seasonal_cycle[mon-1, 3] = np.nanmean(data_temp[qbm_lfm[sel]])

        # if b == 6:
        #     seasonal_cycle[mon-1, b] = np.nanmean(enso.data[sel].copy())
        #     plt.text(11.05, seasonal_cycle[-1, b], lab, horizontalalignment = 'left', verticalalignment = 'baseline', fontsize = 19)
        #     # seasonal_cycle[mon-1, b] = np.nanmean(ac[sel].copy())
        # elif b < 3:
        #     data_temp = enso.data[sel].copy()
        #     # data_temp = ac[sel].copy()
        #     sel2 = (bins_ts[sel] == b)
        #     seasonal_cycle[mon-1, b] = np.nanmean(data_temp[sel2])
        #     plt.text(11.05, seasonal_cycle[-1, b] + 0.07, lab, horizontalalignment = 'left', verticalalignment = 'center', fontsize = 19, 
        #         fontdict = {'color' : col})
        # else:
        #     data_temp = enso.data[sel].copy()
        #     # data_temp = ac[sel].copy()
        #     sel2 = (bins_ts2[sel] == b - 3)
        #     seasonal_cycle[mon-1, b] = np.nanmean(data_temp[sel2])
        #     plt.text(11.05, seasonal_cycle[-1, b] - 0.07, lab, horizontalalignment = 'left', verticalalignment = 'center', fontsize = 19, 
        #         fontdict = {'color' : col})

plt.plot(seasonal_cycle[:, 0], linewidth = 2.)
plt.text(10.05, seasonal_cycle[-1, 0], "QB+\LF+", horizontalalignment = 'left', verticalalignment = 'center', fontsize = 14)
plt.plot(seasonal_cycle[:, 1], linewidth = 2.)
plt.text(10.05, seasonal_cycle[-1, 1], "QB+\LF-", horizontalalignment = 'left', verticalalignment = 'center', fontsize = 14)
plt.plot(seasonal_cycle[:, 2], linewidth = 2.)
plt.text(10.05, seasonal_cycle[-1, 2], "QB-\LF+", horizontalalignment = 'left', verticalalignment = 'center', fontsize = 14)
plt.plot(seasonal_cycle[:, 3], linewidth = 2.)
plt.text(10.05, seasonal_cycle[-1, 3], "QB-\LF-", horizontalalignment = 'left', verticalalignment = 'center', fontsize = 14)
# #     # plt.text(0.8, 0.75, "%.3f$\pm$%.3f" % (np.mean(ac[ndx]), np.std(ac[ndx], ddof = 1)), horizontalalignment = 'center', 
# #         # verticalalignment = 'center', transform = plt.gca().transAxes, fontsize = 14)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xticks(np.arange(0,12), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], size = 15, rotation = 30)
plt.ylabel("NINO3.4 [$^\circ$C]", size = 20)
plt.yticks(size = 15)
plt.ylim([22,29])

ax = plt.gca().twinx()

plt.plot(qbp_lfp_list, ":", color = "C0")
plt.plot(qbp_lfm_list, ":", color = "C1")
plt.plot(qbm_lfp_list, ":", color = "C2")
plt.plot(qbm_lfm_list, ":", color = "C3")
plt.ylim([0,30])
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)

# # plt.subplot(212)
# seasonal = np.zeros((12, 3))
# for mon in range(1,13):
#     ndx = enso.select_months(months = [mon], apply_to_data = False)
#     seasonal[mon-1, 0] = np.nanmean(qbplus[ndx]) # positive
#     seasonal[mon-1, 1] = np.nanmean(qbminus[ndx]) # negative
#     seasonal[mon-1, 2] = np.nanstd(qb[ndx]) # QB variance

# ax.plot(seasonal[:, 0], color = "#c03634", linewidth = 1.6)
# ax.text(11.05, seasonal[-1, 0], "QB+", horizontalalignment = 'left', verticalalignment = 'center', fontsize = 19, 
#                 fontdict = {'color' : "#c03634"})
# ax.plot(seasonal[:, 1], color = "#1f77b4", linewidth = 1.6)
# ax.text(11.05, seasonal[-1, 1], "QB-", horizontalalignment = 'left', verticalalignment = 'center', fontsize = 19, 
#                 fontdict = {'color' : "#1f77b4"})
# ax.plot(seasonal[:, 2], color = "k", linewidth = 1, label = "QB STD")
# ax.text(11.05, seasonal[-1, 2], "QB STD", horizontalalignment = 'left', verticalalignment = 'center', fontsize = 19, 
#                 fontdict = {'color' : "k"})
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.gca().spines['left'].set_visible(False)
# plt.xticks(np.arange(0,12), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], size = 15, rotation = 30)
# plt.ylabel("QB [$^\circ$C]", size = 20)
# plt.yticks(size = 15)
# plt.legend(prop = {'size' : 10})
plt.show()


