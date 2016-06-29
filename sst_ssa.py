import numpy as np
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
from src.data_class import load_ERSST_data
from src.ssa import m_ssa
from datetime import date
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


sst = load_ERSST_data('/Users/nikola/work-ui/data/ersstv4/', date(1900, 1, 1), date(2011,1,1), [-10, 10], [170, 270], True)
if sst.nans:
    sst.interpolate_spatial_nans(method = 'cubic', apply_to_data = True)
sst.subsample_spatial(lat_to = 2, lon_to = 4, start = [0, 170], average = True)
# sst_ts = np.mean(sst.data, axis = 1)
sst.flatten_field()
sst_ts = sst.data.copy()
print sst_ts.shape # averaged over lat -10 - 10, lon by 8 degrees


M = 61
sst_ssa = m_ssa(sst_ts, M = M) # 61 months
lam, e, pc, rc = sst_ssa.run_ssa()
print lam.shape, e.shape, pc.shape, rc.shape
lr, er, pcr, rcr = sst_ssa.apply_varimax(S = 20, sort_lam = True)
print lr.shape, er.shape, pcr.shape, rcr.shape



plt.semilogy(lam[:60], 's', linestyle = 'none', fillstyle = 'none', markersize = 10, color = 'b', label = "SST unrot")
plt.semilogy(lr[:60], linestyle = 'none', marker = 'o', fillstyle = 'none', color = 'k', markersize = 10, label = "SST rot")
plt.legend(loc = 1)
plt.savefig("plots/SSA-SST-eigenvals-hires.png")
plt.show()

# import scipy.signal as ss
# # unrot
# f, px = ss.welch(pc[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, '-', color = '#241632', linewidth = 1.2, label = "PC2 unrot")
# f, px = ss.welch(pc[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, '-', color = '#DE5F48', linewidth = 1.2, label = "PC4 unrot")
# f, px = ss.welch(pc[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, '-', color = '#C7BE2E', linewidth = 1.2, label = "PC6 unrot")

# # rot
# f, px = ss.welch(pcr[:, 1], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, '--', color = '#241632', linewidth = 1.2, label = "PC2 rot")
# f, px = ss.welch(pcr[:, 3], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, '--', color = '#DE5F48', linewidth = 1.2, label = "PC4 rot")
# f, px = ss.welch(pcr[:, 5], fs = 1./2.628e+6, scaling = 'spectrum')
# f *= 3.154e+7
# plt.semilogy(f, px, '--', color = '#C7BE2E', linewidth = 1.2, label = "PC6 rot")

# plt.legend()
# plt.xlim([0, 2])
# plt.gca().xaxis.set_major_locator(MultipleLocator(0.5))
# plt.gca().xaxis.set_minor_locator(MultipleLocator(0.2))

# plt.axvline(1, linestyle = '-.', color = "0.7", linewidth = 0.9)
# plt.axvline(0.5, linestyle = '-.', color = "0.7", linewidth = 0.9)
# plt.axvline(0.25, linestyle = '-.', color = "0.7", linewidth = 0.9)

# plt.ylabel("spectral power")
# plt.xlabel("cycles/year")
# plt.savefig("plots/SSA-SST-PCspectra.png")



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