"""
Main script for ploting CMI maps.
"""
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from main_compute_CMImap import PERIODS_BIN, PERIODS_SPAN, SAVING_FILENAME
from matplotlib.ticker import FuncFormatter, MultipleLocator
from tools import ResultsContainer, SurrogatesContainer

TITLES = ['PHASE SYNCHRONIZATION', 'PHASE-PHASE CAUSALITY', 'PHASE x AMP MI',
          'PHASE-AMP CAUSALITY']
LABELS = ['PHASE', 'PHASE', 'AMP', 'AMP']
ALGORITHM = 'knn'


def get_p_values(data_result, surrogate_result):
    """
    return field with p-values.
    surrogate - (x, y, num_surr)
    """
    assert data_result.shape == surrogate_result.shape[:-1], 'Wrong shapes'
    num_surr = surrogate_result.shape[-1]
    p_vals = np.zeros_like(data_result)
    for i in range(data_result.shape[0]):
        for j in range(data_result.shape[1]):
            p_vals[i, j] = np.sum(np.greater(data_result[i, j],
                                             surrogate_result[i, j, :]))
    p_vals /= float(num_surr)
    return p_vals


def main():
    # load data results
    data_results = ResultsContainer.from_saved_file(
        filename=SAVING_FILENAME + '_data.bin')
    # load surrogate results
    surrogate_result = SurrogatesContainer.from_saved_file(
        filename=SAVING_FILENAME + '_surrogates.bin')

    # plot
    scales = np.arange(PERIODS_SPAN[0], PERIODS_SPAN[-1] + 1, PERIODS_BIN)
    x, y = np.meshgrid(scales, scales)

    plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(2, 2)
    gs.update(left=0.05, right=0.95, hspace=0.3,
              top=0.95, bottom=0.05, wspace=0.15)
    axs = [gs[0, 0], gs[0, 1], gs[1, 0], gs[1, 1]]

    plot = [get_p_values(data_results[result_type, ALGORITHM],
                         surrogate_result[result_type, ALGORITHM]).transpose()
            for result_type in data_results.RESULT_TYPES]

    for ax, cont, tit, lab in zip(axs, plot, TITLES, LABELS):
        ax = plt.subplot(ax)
        ax.contourf(x, y, cont, levels=np.arange(0.95, 1, 0.00125),
                    cmap=plt.cm.get_cmap("jet"), extend='max')
        ax.tick_params(axis='both', which='major', labelsize=17)
        ax.set_title(tit, size=28)
        ax.xaxis.set_major_locator(MultipleLocator(12))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x) / 12))
        ax.xaxis.set_minor_locator(MultipleLocator(6))
        ax.yaxis.set_major_locator(MultipleLocator(12))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x) / 12))
        ax.yaxis.set_minor_locator(MultipleLocator(6))
        ax.set_xlabel("PERIOD PHASE [years]", size=23)
        # plt.colorbar(cs)
        ax.grid()
        ax.set_ylabel("PERIOD %s [years]" % lab, size=23)

    plt.savefig("new_plots/tas_Amon_MPI-ESM-HR_%s.png" % ALGORITHM,
                bbox_inches="tight")


if __name__ == "__main__":
    main()
