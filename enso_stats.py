import numpy as np
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
sys.path.append('/home/nikola/Work/multi-scale')
from src.data_class import load_enso_index, DataField, load_ERSST_data
from datetime import date
import matplotlib.pyplot as plt
import scipy.io as sio

n34 = load_enso_index('nino34raw.txt', 'n34', date(1900,1,1), date(2000,1,1), True)
# n34.get_seasonality()
# sst = load_ERSST_data('/home/nikola/Work/climate-data/ersstv4/', date(1870, 1, 1), date(2006,1,1), [-5, 5], [190, 240], True)
# n34_ersst = np.mean(sst.data, axis=(1,2))

def plot_cmip(enso, models, emr, title = '', ylab = '', fname = None):
    plt.figure(figsize=(12,6))
    for i, model in zip(range(1, len(models)+1), models):
        plt.plot(i, model, 'o', color = '#AB761A', markersize = 10)
    for i, e in zip(range(len(models)+6, len(models)+6 + len(emr)+1), emr):
        plt.plot(i, e, 'o', color = '#9F1C08', markersize = 10)
    plt.axhline(enso, 0, 1, color = '#2F1F21', linewidth = 4, label = 'HadISST N34 index')
    plt.errorbar(0, np.mean(models), yerr = np.std(models, ddof = 1), color = '#AB761A', marker = 's', markersize = 12, label = 'CMIP5 N34')
    plt.errorbar(len(models)+5, np.mean(emr), yerr = np.std(emr, ddof = 1), color = '#9F1C08', marker = 's', markersize = 12, label = 'EMR N34')
    plt.xlim(-1, len(models)+1+len(emr)+5)
    # plt.xlabel('CMIP5 model ensemble mean')
    mds = ['cond. STD 5d %d' % (i+1) for i in range(len(emr))]
    # plt.xticks(np.arange(0, len(models)+1+len(emr)+5), ['mean CMIP'] + model_names + ['', '', '', ''] + ['mean stat.'] + mds, rotation = 75)
    model_names = ['linear %d' % (i+1) for i in range(len(models))]
    plt.xticks(np.arange(0, len(models)+1+len(emr)+5), ['mean linear'] + model_names + ['', '', '', ''] + ['mean cond. STD 5d'] + mds, rotation = 75)
    plt.ylabel(ylab)
    plt.title(title, size = 20)
    # plt.legend(loc=4)
    plt.gca().set_frame_on(False)
    plt.gca().get_xaxis().tick_bottom()
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = 'tight')
    plt.close()


suff = "eps"

model_names = ['CanESM2', 'CCSM4', 'CNRMCM5', 'CSIROmk360', 'GFDLCM3', 'GISSE2Hp1', 'GISSE2Hp2', 'GISSE2Hp3',
    'GISSE2Rp1', 'GISSE2Rp2', 'GISSE2Rp3', 'HadGem2ES', 'IPSL_CM5A_LR', 'MIROC5', 'MRICGCM3']
# models = ['HadGem2ES', 'GISSE2Hp1', 'IPSL_CM5A_LR']
cmip = {}
for model in model_names:
    g = DataField()
    raw = np.loadtxt("N34_CMIP5/N34_%s.txt" % model)
    raw = raw[-1332:, :]
    g.time = n34.time.copy()
    g.data = np.mean(raw, axis = 1)
    g.anomalise()
    cmip[model] = g
    # g.smoothing_running_avg(points = 24, cut_edges = False, use_to_data = True)
    # g.get_seasonality()

emr = {}
# raw = sio.loadmat("Nino34-ERM-1884-2013quadratic-21PCs-sigma0.2.mat")['N34s']
raw = sio.loadmat("S-Nino34-cond-noise-std-5days.mat")['N34s']
for i in range(100):
# for i in range(256,276):
    g = DataField()
    g.time = n34.time.copy()
    g.data = raw[:, i]
    g.anomalise()
    emr[i] = g

emr2 = {}
raw = sio.loadmat("Nino34-ERM-1884-2013linear.mat")['N34s']
for i in range(20):
# for i in range(256,276):
    g = DataField()
    g.time = n34.time.copy()
    g.data = raw[192:-168, i]
    g.anomalise()
    emr2[i] = g

amp = []
ampe = []
# for model in cmip:
#     amp.append(np.std(cmip[model].data, ddof = 1))
# for e in emr:
#     ampe.append(np.std(emr[e].data, ddof = 1))
# for e in emr2:
#     amp.append(np.std(emr2[e].data, ddof = 1))

# plot_cmip(np.std(n34.data, ddof = 1), amp, ampe, 'ENSO amplitude', 'STD of SSTA -- Nino3.4', 'plots/CMIPstats/pub/ENSOampPNF-std-5days.' + suff)

# import scipy.signal as sts
# freq = []
# ps = []
# ps2 = []
# freq2 = []
# f, p = sts.welch(n34.data, 1./2.628e+6)
# f *= 3.154e+7
# ndx1 = np.logical_and(f <= 1/3., f >= 1/8.)
# ndx2 = np.logical_and(f >= 1/3., f <= 1.)
# enso_f = np.mean(p[ndx1]) / np.mean(p[ndx2])
# plt.figure()
# plt.loglog(f, p, linewidth = 3.5, color = '#2F1F21', label = 'HadISST N34 index', zorder = 0)
# plt.xlim(f[0],f[-1])
# plt.xlabel('period [1/yr]')
# plt.ylabel('power')

# for model in cmip:
#     f, p = sts.welch(cmip[model].data, 1./2.628e+6)
#     f *= 3.154e+7
#     freq.append(np.mean(p[ndx1]) / np.mean(p[ndx2]))
#     plt.loglog(f, p, linewidth = 0.5, color = '#9F1C08')
#     ps.append(p)

# for e in emr:
#     f, p = sts.welch(emr[e].data, 1./2.628e+6)
#     f *= 3.154e+7
#     freq2.append(np.mean(p[ndx1]) / np.mean(p[ndx2]))
#     plt.loglog(f, p, linewidth = 0.1, color = '#9F1C08', alpha = 0.6)
#     ps2.append(p)

# for e in emr2:
#     f, p = sts.welch(emr2[e].data, 1./2.628e+6)
#     f *= 3.154e+7
#     freq2.append(np.mean(p[ndx1]) / np.mean(p[ndx2]))
#     plt.loglog(f, p, linewidth = 0.1, color = '#AB761A', alpha = 0.6)
#     ps.append(p)

# ps = np.array(ps)
# ps2 = np.array(ps2)
# plt.loglog(f, np.mean(ps, axis = 0), linewidth = 3, color = '#AB761A', label = 'CMIP5 N34', alpha = 0.6)
# plt.loglog(f, np.mean(ps2, axis = 0), linewidth = 3, color = '#9F1C08', label = 'EMR N34', alpha = 0.6)
# # plt.legend(loc = 3)
# plt.gca().set_frame_on(False)
# plt.gca().get_xaxis().tick_bottom()
# plt.savefig('plots/CMIPstats/pub/ENSOspectraPNF-std-5days.' + suff, bbox_inches = 'tight')
# plt.close()

# plot_cmip(enso_f, freq, freq2, 'ENSO frequency', 'P(3-8years) / P(1-3years)', 'plots/CMIPstats/pub/ENSOfreqPNF-std-5days.' + suff)

mons = list(np.arange(1,13))
enso_s = []
for mon in mons:
    ndx = n34.select_months([mon], False)
    enso_s.append(np.std(n34.data[ndx], ddof = 1))
plt.plot(mons, enso_s, linewidth = 3.5, color = '#2F1F21', label = 'HadISST N34 index')
plt.xlim(1,12)
plt.xlabel('month')
plt.ylabel('STD of SSTA')
enso_s = np.array(enso_s)
ndj_ndx = [0,10,11]
mam_ndx = [2,3,4]
seas = []
model_s = []
seas2 = []
model_s2 = []
# for model in cmip:
#     mod = []
#     for mon in mons:
#         ndx = n34.select_months([mon], False)
#         mod.append(np.std(cmip[model].data[ndx], ddof = 1))
#     plt.plot(mons, mod, linewidth = 0.5, color = '#9F1C08')
#     seas.append(mod)
#     mod = np.array(mod)
#     model_s.append(np.mean(mod[ndj_ndx]) / np.mean(mod[mam_ndx]))

for e in emr:
    mod = []
    for mon in mons:
        ndx = n34.select_months([mon], False)
        mod.append(np.std(emr[e].data[ndx], ddof = 1))
    plt.plot(mons, mod, linewidth = 0.5, color = '#9F1C08')
    seas2.append(mod)
    mod = np.array(mod)
    model_s2.append(np.mean(mod[ndj_ndx]) / np.mean(mod[mam_ndx]))

for e in emr2:
    mod = []
    for mon in mons:
        ndx = n34.select_months([mon], False)
        mod.append(np.std(emr2[e].data[ndx], ddof = 1))
    plt.plot(mons, mod, linewidth = 0.5, color = '#AB761A')
    seas.append(mod)
    mod = np.array(mod)
    model_s.append(np.mean(mod[ndj_ndx]) / np.mean(mod[mam_ndx]))

seas = np.array(seas)
seas2 = np.array(seas2)
plt.plot(mons, np.mean(seas, axis = 0), linewidth = 3, color = '#AB761A', label = 'CMIP5 N34')
plt.plot(mons, np.mean(seas2, axis = 0), linewidth = 3, color = '#9F1C08', label = 'EMR N34')
mo = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
plt.xticks(np.arange(1,13), mo, rotation = 45)
plt.gca().set_frame_on(False)
plt.gca().get_xaxis().tick_bottom()
plt.savefig('plots/CMIPstats/pub/ENSOseasonalityPNF-std-5days.' + suff, bbox_inches = 'tight')

# plot_cmip(np.mean(enso_s[ndj_ndx]) / np.mean(enso_s[mam_ndx]), model_s, model_s2, 'ENSO seasonality', 'STD of SSTA (NDJ / MAM)', 'plots/CMIPstats/pub/ENSOseasPNF-std-5days.' + suff)