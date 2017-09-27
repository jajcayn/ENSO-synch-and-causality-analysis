"""
Plots figure with CMIP5 comparison.
"""


import csv
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

with open("CMIP5-interactions.csv", 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    model = []
    datal1 = []
    datal2 = []
    datacorr = []
    data_int = []
    for row in reader:
        if row[0] == 'model':
            continue
        if 'mean' in row[0]:
            model.append(row[0][:-9])
            datal1.append(float(row[1]))
            datal2.append(float(row[2]))
            datacorr.append(float(row[3]))
            data_int.append([float(i) for i in row[4:]])

datal1 = np.array(datal1)[:, np.newaxis]
print datal1.shape
datal2 = np.array(datal2)[:, np.newaxis]
print datal2.shape
datacorr = np.array(datacorr)[:, np.newaxis]
print datacorr.shape
data_int = np.array(data_int)
print data_int.shape


plt.figure()
# gs = gridspec.GridSpec(1, 4, width_ratios = [1, 1, 1, 6])
gs = gridspec.GridSpec(1, 3, width_ratios = [1, 1, 6])
gs.update(wspace = 0.2)
# ax = plt.subplot(gs[0])
# cs = plt.pcolor(datal1, cmap = plt.get_cmap('hot_r'), vmin = 180, vmax = 800)
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom='off',      # ticks along the bottom edge are off
#     top='off',         # ticks along the top edge are off
#     labelbottom='off',
#     labeltop='on')
# plt.tick_params(
#     axis='y',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     right='off',      # ticks along the bottom edge are off,         # ticks along the top edge are off
#     )
# ax.set_xticks([0.5])
# ax.set_xticklabels(['L1 distance'], rotation = 90, size = 20)
# ax.set_aspect('equal')
# ax.set_yticks(np.arange(0.5, 15.5, 1))
# ax.set_yticklabels(model, size = 20)
# ax.set_xlim([0,1])
# ax.set_ylim([0,15])
# plt.colorbar(cs, fraction = 0.15, pad = 0.08, ticks = np.arange(200,1000,200))
ax = plt.subplot(gs[0])
cs = plt.pcolor(datal2, cmap = plt.get_cmap('hot_r'), vmin = 25, vmax = 120)
ax.set_aspect('equal')
ax.set_xlim([0,1])
ax.set_ylim([0,15])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off',
    labeltop='on')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    right='off',      # ticks along the bottom edge are off,         # ticks along the top edge are off
    )
ax.set_yticks(np.arange(0.5, 15.5, 1))
ax.set_yticklabels(model, size = 20)
ax.set_xticks([0.5])
ax.set_xticks([0.5])
ax.set_xticklabels(['L2 distance'], rotation = 90, size = 20)
plt.colorbar(cs, fraction = 0.15, pad = 0.08, ticks = np.arange(30,150,30))
ax = plt.subplot(gs[1])
cs = plt.pcolor(datacorr, cmap = plt.get_cmap('hot'), vmin = 0, vmax = 0.8)
ax.set_aspect('equal')
ax.set_xlim([0,1])
ax.set_ylim([0,15])
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off',
    labeltop='on')
ax.set_xticks([0.5])
ax.set_xticklabels(['correlation'], rotation = 90, size = 20)
ax.set_yticks([], [])
plt.colorbar(cs, fraction = 0.15, pad = 0.08, ticks = np.arange(0,1,0.2))
ax = plt.subplot(gs[2])
cs = plt.pcolor(data_int, cmap = plt.get_cmap('hot'), vmin = 0, vmax = 0.25)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off',
    labeltop='on')
ax.set_xlim([0,6])
ax.set_ylim([0,15])
ax.set_aspect('equal')
ax.set_yticks([], [])
ax.set_xticks(np.arange(0.5,6.5))
ax.set_xticklabels(['SYNCH - corr', 'SYNCH - ARI', 'PHcausality - corr', 'PHcausality - ARI', 'P-Acausality - corr', 'P-Acausality - ARI'], rotation = 90, size = 20)
plt.colorbar(cs, fraction = 0.15*0.15, pad = 0.08, ticks = np.arange(0, 0.3, 0.05))

# plt.show()
plt.savefig('plots/CMIP5_model_ens_test.eps', bbox_inches = 'tight')
