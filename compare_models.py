import numpy as np
import cPickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FuncFormatter
import os

BINS = 4
WVLT_SPAN = [5,93]

def evaluate_MI(fname, threshold = 0.95):
    CUT = slice(0,1000)
    # version = 3
    with open(fname, 'rb') as f:
        result = cPickle.load(f)

    phase_phase_coherence = result['phase x phase data']
    phase_phase_CMI = result['phase CMI data']
    surrCoherence = result['phase x phase surrs'][CUT, ...]
    surrCMI = result['phase CMI surrs'][CUT, ...]
    phase_amp_MI = result['phase x amp data']
    phase_amp_condMI = result['phase amp CMI data']
    surrPhaseAmp = result['phase x amp surrs'][CUT, ...]
    surrPhaseAmpCMI = result['phase amp CMI surrs'][CUT, ...]

    res_phase_coh = np.zeros_like(phase_phase_coherence)
    res_phase_cmi = np.zeros_like(res_phase_coh)
    res_phase_amp = np.zeros_like(res_phase_coh)
    res_phase_amp_CMI = np.zeros_like(res_phase_coh)

    for i in range(res_phase_coh.shape[0]):
        for j in range(res_phase_coh.shape[1]):
            res_phase_coh[i, j] = np.sum(np.greater(phase_phase_coherence[i, j], surrCoherence[:, i, j])) / np.float(surrCoherence.shape[0])
            res_phase_cmi[i, j] = np.sum(np.greater(phase_phase_CMI[i, j], surrCMI[:, i, j])) / np.float(surrCMI.shape[0])
            res_phase_amp[i, j] = np.sum(np.greater(phase_amp_MI[i, j], surrPhaseAmp[:, i, j])) / np.float(surrPhaseAmp.shape[0])
            res_phase_amp_CMI[i, j] = np.sum(np.greater(phase_amp_condMI[i, j], surrPhaseAmpCMI[:, i, j])) / np.float(surrPhaseAmpCMI.shape[0])

    f.close()

    res_phase_coh_thr = np.zeros_like(res_phase_coh, dtype = np.int)
    res_phase_coh_thr[np.where(res_phase_coh > threshold)] = 1
    res_phase_cmi_thr = np.zeros_like(res_phase_cmi, dtype = np.int)
    res_phase_cmi_thr[np.where(res_phase_cmi > threshold)] = 1
    res_phase_amp_thr = np.zeros_like(res_phase_amp, dtype = np.int)
    res_phase_amp_thr[np.where(res_phase_amp > threshold)] = 1
    res_phase_amp_CMI_thr = np.zeros_like(res_phase_amp_CMI, dtype = np.int)
    res_phase_amp_CMI_thr[np.where(res_phase_amp_CMI > threshold)] = 1

    return res_phase_coh_thr, res_phase_cmi_thr, res_phase_amp_thr, res_phase_amp_CMI_thr


def plot_MI(MIs, data, fname = "test.png"):
    scales = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)
    # a = np.random.rand(scales.shape[0], scales.shape[0]) + 0.5
    x, y = np.meshgrid(scales, scales)

    # fig, axs = plt.subplots(1, 2, figsize = (13,7))
    fig = plt.figure(figsize=(15,15))
    gs = gridspec.GridSpec(2, 2)
    gs.update(left=0.05, right=0.95, hspace=0.3, top=0.95, bottom=0.05, wspace=0.15)
    i = 0
    axs = [gs[0,0], gs[0,1], gs[1,0], gs[1,1]]
    # plot = [res_phase_coh.T, res_phase_cmi.T, res_phase_amp.T, res_phase_amp_CMI.T]
    plot = [MIs[0].T, MIs[1].T, MIs[2].T, MIs[3].T]
    tits = ['PHASE COHERENCE', 'CMI PHASE DIFF', 'PHASE x AMP MI', 'PHASE x AMP CMI Gauss']
    counts = []
    for ax, cont, tit, dat in zip(axs, plot, tits, data):
        ax = plt.subplot(ax)
        cs = ax.contourf(x, y, cont, 2, cmap = plt.cm.get_cmap("jet"), extend = 'max')
        ax.tick_params(axis='both', which='major', labelsize = 17)
        ax.set_title(tit, size = 28)
        test_string = ("%d // %d -- miss %d" % (cont[cont == 2].shape[0], dat[dat == 1].shape[0], cont[cont == 1].shape[0]))
        counts.append([cont[cont == 2].shape[0], dat[dat == 1].shape[0], cont[cont == 1].shape[0]])
        ax.text(0.62, 0.9, test_string, horizontalalignment = 'center', verticalalignment = 'center', size = 27, 
            transform = ax.transAxes, color = "white", bbox = dict(boxstyle = "square,pad=0.25", facecolor = 'grey', alpha = 0.7))
        ax.xaxis.set_major_locator(MultipleLocator(12))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
        ax.xaxis.set_minor_locator(MultipleLocator(6))
        ax.yaxis.set_major_locator(MultipleLocator(12))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
        ax.yaxis.set_minor_locator(MultipleLocator(6))
        ax.set_xlabel("period [years]", size = 20)
        if i % 2 == 0:
            ax.set_ylabel("period [years]", size = 20)
        else:
            # fig.colorbar(cs, ax = ax, shrink = 0.5)
            pass
        i += 1

    plt.savefig(fname)
    plt.close()

    return counts


## data
fname = ("CMImap%dbins3Dcond_GaussCorr.bin" % (BINS))
data = evaluate_MI(fname)

# ## models
# for a in os.walk("models"):
#     if ".bin" in a[2][0]:
#         break

# models = a[2]

CMIP5models = ['N34_CanESM2', 'N34_GFDLCM3', 'N34_GISSE2Hp1', 'N34_GISSE2Hp2', 'N34_GISSE2Hp3', 'N34_GISSE2Rp1']
CMIP5models += ['N34_GISSE2Rp2', 'N34_GISSE2Rp3', 'N34_HadGem2ES', 'N34_IPSL_CM5A_LR', 'N34_MIROC5', 'N34_MRICGCM3']
CMIP5models += ['N34_CCSM4', 'N34_CNRMCM5', 'N34_CSIROmk360']

numbers = []

for CMIP5model in CMIP5models:
    fname = CMIP5model + '.txt'
    model = np.loadtxt('N34_CMIP5/' + fname)
    model_count = model.shape[1]

    for num_ts in range(model_count):
        result_temp = evaluate_MI("models/CMImap%dbins3Dcond_GaussCorr_%sts%d.bin" % (BINS, CMIP5model, num_ts))
        print("%s time series %d model read.." % (CMIP5model, num_ts+1))

        # tests
        test = []
        for no in range(4):
            summed = result_temp[no] + data[no]
            test.append(summed)

        fname = ("test%sts%d.png" % (CMIP5model, num_ts))
        counts = plot_MI(test, data, fname = "models/plots/" + fname)
        numbers.append(counts)

numbers = np.array(numbers)
print numbers.shape

# idx = numbers[:, 0, 0].argmax()
# print("Highest agreement PHASE COHERENCE: %s with %d // %d -- %.1f%%" % (CMIP5models[idx], 
#     numbers[idx, 0, 0], numbers[idx, 0, 1], numbers[idx, 0, 0] / float(numbers[idx, 0, 1]) * 100))

# idx = numbers[:, 1, 0].argmax()
# print("Highest agreement CMI PHASE DIFF: %s with %d // %d -- %.1f%%" % (CMIP5models[idx], 
#     numbers[idx, 1, 0], numbers[idx, 1, 1], numbers[idx, 1, 0] / float(numbers[idx, 1, 1]) * 100))

# idx = numbers[:, 2, 0].argmax()
# print("Highest agreement PHASE x AMP MI: %s with %d // %d -- %.1f%%" % (CMIP5models[idx], 
#     numbers[idx, 2, 0], numbers[idx, 2, 1], numbers[idx, 2, 0] / float(numbers[idx, 2, 1]) * 100))

# idx = numbers[:, 3, 0].argmax()
# print("Highest agreement PHASE x AMP CMI Gauss: %s with %d // %d -- %.1f%%" % (CMIP5models[idx], 
#     numbers[idx, 3, 0], numbers[idx, 3, 1], numbers[idx, 3, 0] / float(numbers[idx, 3, 1]) * 100))
