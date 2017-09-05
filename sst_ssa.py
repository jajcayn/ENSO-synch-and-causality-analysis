import numpy as np
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
from src.data_class import load_ERSST_data, load_enso_index
from src.ssa import ssa_class
from datetime import date
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


sst = load_ERSST_data('/Users/nikola/work-ui/data/ersstv4/', date(1900, 1, 1), date(2011,1,1), [-8, 8], [170, 270], True)
if sst.nans:
    sst.interpolate_spatial_nans(method = 'cubic', apply_to_data = True)
sst.subsample_spatial(lat_to = 2, lon_to = 4, start = [0, 170], average = True)
# sst_ts = np.mean(sst.data, axis = 1)
print sst.data.shape
sst.flatten_field()
sst_ts = sst.data.copy()
print sst_ts.shape # averaged over lat -10 - 10, lon by 8 degrees

nino34 = load_enso_index("/Users/nikola/work-ui/data/nino34raw.txt", '3.4', date(1900, 1, 1), date(2011, 1, 1), anom = True)


M = 72
sst_ssa = ssa_class(sst_ts, M = M, compute_rc = False) # 61 months
nino_ssa = ssa_class(nino34.data, M = M, compute_rc = False)
lam, e, pc = sst_ssa.run_ssa()
lam_nino, e_nino, pc_nino = nino_ssa.run_ssa()
# num = sst_ssa.run_Monte_Carlo(1000, method = "rank-def")
# print lam.shape, e.shape, pc.shape
lr, er, pcr = sst_ssa.apply_varimax(S = 20, sort_lam = True)
# lr_nino, er_nino, pcr_nino = nino_ssa.apply_varimax(S = 20, sort_lam = True)
print lr.shape, er.shape, pcr.shape


plt.semilogy(lam[:40], 's', linestyle = 'none', fillstyle = 'none', markersize = 10, color = 'b', label = "SST unrot")
plt.semilogy(lr[:40], linestyle = 'none', marker = 'o', fillstyle = 'none', color = 'k', markersize = 10, label = "SST rot")
plt.semilogy(lam_nino[:40], 'v', linestyle = 'none', fillstyle = 'none', markersize = 10, color = 'b', label = "NINO")
# plt.semilogy(lr_nino[:40], linestyle = 'none', marker = '^', fillstyle = 'none', color = 'k', markersize = 10, label = "NINO rot")
plt.legend(loc = 1)
# plt.savefig("plots/SSA-SST-eigenvals-hires.png")
plt.show()

import scipy.signal as ss
plt.figure()
# # unrot
f, px = ss.welch(pc[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-', color = '#241632', linewidth = 1.2, label = "SST: PC2 unrot")
f, px = ss.welch(pc[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-', color = '#DE5F48', linewidth = 1.2, label = "SST: PC4 unrot")
f, px = ss.welch(pc[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-', color = '#C7BE2E', linewidth = 1.2, label = "SST: PC6 unrot")
#
f, px = ss.welch(pc_nino[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-.', color = '#241632', linewidth = 1.2, label = "NINO: PC2")
f, px = ss.welch(pc_nino[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-.', color = '#DE5F48', linewidth = 1.2, label = "NINO: PC4")
f, px = ss.welch(pc_nino[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-.', color = '#C7BE2E', linewidth = 1.2, label = "NINO: PC6")

# rot
f, px = ss.welch(pcr[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '--', color = '#241632', linewidth = 1.2, label = "SST: PC2 rot")
f, px = ss.welch(pcr[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '--', color = '#DE5F48', linewidth = 1.2, label = "SST: PC4 rot")
f, px = ss.welch(pcr[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '--', color = '#C7BE2E', linewidth = 1.2, label = "SST: PC6 rot")
#
# f, px = ss.welch(pcr_nino[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, ':', color = '#241632', linewidth = 1.2, label = "NINO: PC2 rot")
# f, px = ss.welch(pcr_nino[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, ':', color = '#DE5F48', linewidth = 1.2, label = "NINO: PC4 rot")
# f, px = ss.welch(pcr_nino[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, ':', color = '#C7BE2E', linewidth = 1.2, label = "NINO: PC6 rot")

plt.legend()
plt.xlim([0, 2])
plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
plt.gca().xaxis.set_minor_locator(MultipleLocator(0.2))

plt.axvline(1, linestyle = '-.', color = "0.7", linewidth = 0.9)
plt.axvline(0.5, linestyle = '-.', color = "0.7", linewidth = 0.9)
plt.axvline(0.25, linestyle = '-.', color = "0.7", linewidth = 0.9)

plt.ylabel("spectral power")
plt.xlabel("cycles/year")
# plt.savefig("plots/SSA-SST-PCspectra246-hires.png")
plt.show()
# plt.close()


plt.figure()
# unrot
f, px = ss.welch(pc[:, 0], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-', color = '#241632', linewidth = 1.2, label = "SST: PC1 unrot")
f, px = ss.welch(pc[:, 2], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-', color = '#DE5F48', linewidth = 1.2, label = "SST: PC3 unrot")
f, px = ss.welch(pc[:, 4], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-', color = '#C7BE2E', linewidth = 1.2, label = "SST: PC5 unrot")

f, px = ss.welch(pc_nino[:, 0], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-.', color = '#241632', linewidth = 1.2, label = "NINO: PC1")
f, px = ss.welch(pc_nino[:, 2], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-.', color = '#DE5F48', linewidth = 1.2, label = "NINO: PC3")
f, px = ss.welch(pc_nino[:, 4], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '-.', color = '#C7BE2E', linewidth = 1.2, label = "NINO: PC5")

# rot
f, px = ss.welch(pcr[:, 0], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '--', color = '#241632', linewidth = 1.2, label = "SST: PC1 rot")
f, px = ss.welch(pcr[:, 2], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '--', color = '#DE5F48', linewidth = 1.2, label = "SST: PC3 rot")
f, px = ss.welch(pcr[:, 4], fs = 1./2.628e+6, scaling = 'spectrum')
f *= 3.154e+7
plt.semilogy(f, px, '--', color = '#C7BE2E', linewidth = 1.2, label = "SST: PC5 rot")

plt.legend()
plt.xlim([0, 2])
plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
plt.gca().xaxis.set_minor_locator(MultipleLocator(0.2))

plt.axvline(1, linestyle = '-.', color = "0.7", linewidth = 0.9)
plt.axvline(0.5, linestyle = '-.', color = "0.7", linewidth = 0.9)
plt.axvline(0.25, linestyle = '-.', color = "0.7", linewidth = 0.9)

plt.ylabel("spectral power")
plt.xlabel("cycles/year")
plt.show()
# plt.savefig("plots/SSA-SST-PCspectra135-hires.png")


# import scipy.io as sio
# e = sio.loadmat('SSA-SST-hires.mat')['e']
# er = sio.loadmat('SSA-SST-hires.mat')['er']

# print e.shape, er.shape

# t = e[:, 2].reshape((M, -1), order = 'F')
# print t.shape
# print sst.data.shape
# t = sst.reshape_flat_field(f = t)
# print t.shape
# # print t.min(), t.max()


# def bring_to(f, to):
#     m = (f.max() - f.min()) / 2.
#     f -= (m - to)

#     return f

# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# x, y = np.meshgrid(sst.lats, sst.lons, indexing = 'ij')
# for i in range(t.shape[0]):
#     ax.contourf(x, y, bring_to(t[i, ...], i), zdir = 'z', offset = t[i, ...].max(), alpha = 0.5)
# ax.set_ylabel("LON")
# ax.set_xlabel("LAT")
# plt.show()

# f, axarr = plt.subplots(1, 4, sharey = True)
# for i in range(4):
#     axarr[i].contourf(e[:, i].reshape((M, -1), order = 'F'), 15, cmap = plt.get_cmap('Greys'))
#     axarr[i].set_xticks(np.arange(0, sst.lons.shape[0],4))
#     axarr[i].set_xticklabels(sst.lons[::4])
#     axarr[i].set_yticks(np.arange(0,M,3))
#     # axarr[i].set_yticklabels(np.arange(0,M,3))
#     if i%2 == 0:
#         axarr[i].set_ylabel("lag [months]")
#     axarr[i].set_xlabel("lon")
#     axarr[i].set_title("SST $\mathbf{e}_{%d}$" % (i+1))

# plt.suptitle("SST SSA eigenvectors: unrotated")

# plt.savefig("plots/SSA-SST-eigenvec-short.png")