import numpy as np
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
from src.data_class import load_ERSST_data, load_enso_index, DataField
from src.ssa import ssa_class
from datetime import date
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import scipy.signal as ss
from parametric_recharge_oscillator import ENSO_PROmodel


CMIP5models = ['N34_CanESM2', 'N34_GFDLCM3', 'N34_GISSE2Hp1', 'N34_GISSE2Hp2', 'N34_GISSE2Hp3', 'N34_GISSE2Rp1']
CMIP5models += ['N34_GISSE2Rp2', 'N34_GISSE2Rp3', 'N34_HadGem2ES', 'N34_IPSL_CM5A_LR', 'N34_MIROC5', 'N34_MRICGCM3']
CMIP5models += ['N34_CCSM4', 'N34_CNRMCM5', 'N34_CSIROmk360']


# sst = load_ERSST_data('/Users/nikola/work-ui/data/ersstv4/', date(1900, 1, 1), date(2011,1,1), [-8, 8], [170, 270], True)
# if sst.nans:
#     sst.interpolate_spatial_nans(method = 'cubic', apply_to_data = True)
# sst.subsample_spatial(lat_to = 2, lon_to = 4, start = [0, 170], average = True)
# # sst_ts = np.mean(sst.data, axis = 1)
# print sst.data.shape
# sst.flatten_field()
# sst_ts = sst.data.copy()
# print sst_ts.shape # averaged over lat -10 - 10, lon by 8 degrees

pro_neutral = DataField()
nino34 = load_enso_index("/Users/nikola/work-ui/data/nino34raw.txt", '3.4', date(1900, 1, 1), date(2005, 1, 1), anom = True)
# PROmodel_enso = ENSO_PROmodel(length = nino34.data.shape[0], daily = False, damped = False, ENSOperiod = 3.75, modulation = 2., lambda0 = 0.4)
# PROmodel_enso.integrate_PROmodel()
# pro_neutral.data = PROmodel_enso.data.copy()
# pro_neutral.create_time_array(date_from = date(1900, 1, 1), sampling = 'm')


M = 72
cmips_eigs = []
cmips_pcs = []
raw = np.loadtxt("Sergey-NINO34-SSTlinear-20PCs-L3-noise.txt")
print raw.shape
for i in range(raw.shape[1]):
    pro_damped = DataField()
    # PROmodel_enso = ENSO_PROmodel(length = nino34.data.shape[0], daily = False, damped = True, ENSOperiod = 3.75, modulation = 2., lambda0 = 0.4)
    # PROmodel_enso.integrate_PROmodel()
    # pro_damped.data = PROmodel_enso.data.copy()
    pro_damped.data = raw[:, i]
    pro_damped.create_time_array(date_from = date(1900, 1, 1), sampling = 'm')
    pro_damped.anomalise()

    pro_ssa = ssa_class(pro_damped.data, M = M, compute_rc = False)
    lam, e, pc = pro_ssa.run_ssa()
    cmips_eigs.append(lam[:40])
    cmips_pcs.append(pc[:, :6])

    # fname = CMIP5model + '.txt'
    # model = np.loadtxt('N34_CMIP5/' + fname)
    # model_count = model.shape[1]

    # print model.shape
    # model_means_eigs = []
    # model_means_pcs = []
    # for num_ts in range(model_count):
#         cmip_ts = DataField(data = model[:, num_ts])
#         cmip_ts.create_time_array(date_from = date(1860, 1, 1), sampling = 'm')
#         cmip_ts.select_date(date(1900, 1, 1), date(2005, 1, 1))
#         cmip_ts.anomalise()
#         cmip_ssa = ssa_class(cmip_ts.data, M = M, compute_rc = False)
#         lam, e, pc = cmip_ssa.run_ssa()
#         model_means_eigs.append(lam[:40])
#         model_means_pcs.append(pc[:, :6])
#     model_means_eigs = np.array(model_means_eigs)
#     model_means_pcs = np.array(model_means_pcs)
#     cmips_eigs.append(np.mean(model_means_eigs, axis = 0))
#     cmips_pcs.append(np.mean(model_means_pcs, axis = 0))
cmips_eigs = np.array(cmips_eigs)
cmips_pcs = np.array(cmips_pcs)

print cmips_eigs.shape, cmips_pcs.shape


# sst_ssa = ssa_class(sst_ts, M = M, compute_rc = False) # 61 months
nino_ssa = ssa_class(nino34.data, M = M, compute_rc = False)
# pro_neutral_ssa = ssa_class(pro_neutral.data, M = M, compute_rc = False)
# lam, e, pc = sst_ssa.run_ssa()
lam_nino, e_nino, pc_nino = nino_ssa.run_ssa()
# lam_pro, e_pro, pc_pro = pro_neutral_ssa.run_ssa()
# lr, er, pcr = sst_ssa.apply_varimax(S = 20, sort_lam = True)


norm = plt.Normalize()
colors = plt.cm.jet(norm(np.arange(raw.shape[1])))

# eigenvalues
# plt.semilogy(lam[:40], 's', linestyle = 'none', fillstyle = 'none', markersize = 10, color = 'C0', label = "SST unrot")
# plt.semilogy(lr[:40], linestyle = 'none', marker = 'o', fillstyle = 'none', color = 'C0', markersize = 10, label = "SST rot")
plt.semilogy(lam_nino[:40], 'v', linestyle = 'none', fillstyle = 'none', markersize = 12, color = 'k', label = "NINO3.4")
# plt.semilogy(lam_pro[:40], 'x', linestyle = 'none', fillstyle = 'none', markersize = 12, color = 'C0', label = "PRO neutral")
for i in range(raw.shape[1]):
    plt.semilogy(cmips_eigs[i, :], '.', linestyle = 'none', fillstyle = 'none', markersize = 5, color = colors[i], 
        label = "empirical")

# plt.legend(loc = 1)
plt.show()

# Welch spectra of PCs
for pc_num in range(6):
    plt.figure()
    # f, px = ss.welch(pc[:, 0], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-', color = 'C0', linewidth = 0.8, label = "SST: PC1 unrot")
    # f, px = ss.welch(pc[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-', color = 'C1', linewidth = 0.8, label = "SST: PC2 unrot")
    # f, px = ss.welch(pc[:, 2], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-', color = 'C2', linewidth = 0.8, label = "SST: PC3 unrot")
    # f, px = ss.welch(pc[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-', color = 'C3', linewidth = 0.8, label = "SST: PC4 unrot")
    # f, px = ss.welch(pc[:, 4], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-', color = 'C4', linewidth = 0.8, label = "SST: PC5 unrot")
    # f, px = ss.welch(pc[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-', color = 'C5', linewidth = 0.8, label = "SST: PC6 unrot")

    f, px = ss.welch(pc_nino[:, pc_num], fs = 1./2.628e+6, scaling = 'spectrum')
    f *= 3.154e+7
    plt.semilogy(f, px, '-', color = 'k', linewidth = 2.5, label = "NINO")

    # f, px = ss.welch(pc_pro[:, pc_num], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '--', color = 'k', linewidth = 2, label = "PRO neutral")
    # f, px = ss.welch(pc_nino[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '--', color = 'k', linewidth = 0.8)#, label = "NINO: PC2")
    # f, px = ss.welch(pc_nino[:, 2], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-.', color = 'k', linewidth = 0.8)#, label = "NINO: PC3")
    # f, px = ss.welch(pc_nino[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, ':', color = 'k', linewidth = 0.8)#, label = "NINO: PC4")
    # f, px = ss.welch(pc_nino[:, 4], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '--', dashes = (5,10), color = 'k', linewidth = 0.8)#, label = "NINO: PC5")
    # f, px = ss.welch(pc_nino[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '--', dashes = (5,20), color = 'k', linewidth = 0.8)#, label = "NINO: PC6")

    for i in range(raw.shape[1]):
        f, px = ss.welch(cmips_pcs[i, :, pc_num], fs = 1./2.628e+6, scaling = 'spectrum')
        f *= 3.154e+7
        plt.semilogy(f, px, '-', color = colors[i], linewidth = 0.5, label = "empirical")
        # f, px = ss.welch(cmips_pcs[i, :, 1], fs = 1./2.628e+6, scaling = 'spectrum')
        # f *= 3.154e+7
        # plt.semilogy(f, px, '--', color = colors[i], linewidth = 0.5)#, label = "%s: PC2" % (CMIP5models[i][4:]))
        # f, px = ss.welch(cmips_pcs[i, :, 2], fs = 1./2.628e+6, scaling = 'spectrum')
        # f *= 3.154e+7
        # plt.semilogy(f, px, '-.', color = colors[i], linewidth = 0.5)#, label = "%s: PC3" % (CMIP5models[i][4:]))
        # f, px = ss.welch(cmips_pcs[i, :, 3], fs = 1./2.628e+6, scaling = 'spectrum')
        # f *= 3.154e+7
        # plt.semilogy(f, px, ':', color = colors[i], linewidth = 0.5)#, label = "%s: PC4" % (CMIP5models[i][4:]))
        # f, px = ss.welch(cmips_pcs[i, :, 4], fs = 1./2.628e+6, scaling = 'spectrum')
        # f *= 3.154e+7
        # plt.semilogy(f, px, '--', dashes = (5,10), color = colors[i], linewidth = 0.5)#, label = "%s: PC5" % (CMIP5models[i][4:]))
        # f, px = ss.welch(cmips_pcs[i, :, 5], fs = 1./2.628e+6, scaling = 'spectrum')
        # f *= 3.154e+7
        # plt.semilogy(f, px, '--', dashes = (5,20), color = colors[i], linewidth = 0.5)#, label = "%s: PC6" % (CMIP5models[i][4:]))

    # f, px = ss.welch(pcr[:, 0], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-.', color = 'C0', linewidth = 0.8, label = "SST: PC1 rot")
    # f, px = ss.welch(pcr[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-.', color = 'C1', linewidth = 0.8, label = "SST: PC2 rot")
    # f, px = ss.welch(pcr[:, 2], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-.', color = 'C2', linewidth = 0.8, label = "SST: PC3 rot")
    # f, px = ss.welch(pcr[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-.', color = 'C3', linewidth = 0.8, label = "SST: PC4 rot")
    # f, px = ss.welch(pcr[:, 4], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-.', color = 'C4', linewidth = 0.8, label = "SST: PC5 rot")
    # f, px = ss.welch(pcr[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
    # f *= 3.154e+7
    # plt.semilogy(f, px, '-.', color = 'C5', linewidth = 0.8, label = "SST: PC6 rot")

    # plt.legend()
    plt.xlim([0, 2])
    plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
    plt.gca().xaxis.set_minor_locator(MultipleLocator(0.2))

    plt.axvline(1, linestyle = '-.', color = "0.7", linewidth = 0.9)
    plt.axvline(0.5, linestyle = '-.', color = "0.7", linewidth = 0.9)
    plt.axvline(0.25, linestyle = '-.', color = "0.7", linewidth = 0.9)

    plt.title("PC%d" % (int(pc_num)+1), size = 25)
    plt.ylabel("spectral power")
    plt.xlabel("cycles/year")
    plt.show()