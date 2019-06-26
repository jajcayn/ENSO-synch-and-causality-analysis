"""
Main script for computing CMI on single timeseries.
"""
from datetime import date, datetime
from multiprocessing import Process, Queue

import numpy as np
import pyclits.mutual_inf as mutual
import xarray as xr
from pyclits.geofield import DataField
from pyclits.surrogates import SurrogateField
from tools import ResultsContainer, SurrogatesContainer
from tqdm import tqdm

PERIODS_SPAN = [6, 240]  # unit is month, 96
PERIODS_BIN = 6  # in months
NUM_SURROGATES = 100
WORKERS = 20
NUM_BINS_EQQ = 4
K_KNN = 64

NINO_INDEX_BOUNDS = {
    "NINO3.4": {"lat": slice(-5, 5), "lon": slice(190, 240)},
    "NINO3": {"lat": slice(-5, 5), "lon": slice(210, 270)},
    "NINO4": {"lat": slice(-5, 5), "lon": slice(160, 210)},
    "NINO1.2": {"lat": slice(-10, 0), "lon": slice(270, 280)}
}

# TO BE CHANGED
DATASET_PATH = 'new_data_may19/nino34.txt'
DATASET_VARIABLE = 'tas'
FIRST_DATE = '1870-01-01'

SAVING_FILENAME = 'bins/nino34_6-240months'


def prepare_dataset(dataset_path, nino_region=None, surrogates=True):
    """
    Prepare dataset into single timeseries.

    :dataset_path: path to the dataset, to be loaded with xarray
    :nino_region: which nino region to use, if None, will not select data
    :surrogates: whether surrogates will be computed
    """
    if dataset_path.endswith('.nc'):
        dataset = xr.open_dataset(dataset_path)
        variable = dataset[DATASET_VARIABLE]
        assert variable.ndim == 3, (
            'Currently supports only geospatial lat-lon datasets')
        # cut NINO region
        if nino_region is not None:
            variable = variable.sel(**NINO_INDEX_BOUNDS[nino_region])
        variable.sel(time=slice(FIRST_DATE, None))
        # do spatial mean
        variable = variable.mean(dim=["lat", "lon"])
        assert variable.ndim == 1, "Now the dataset should be 1dimensional"
        timeseries = DataField(data=variable.values)

    elif dataset_path.endswith('.txt'):
        dataset = np.loadtxt(dataset_path, skiprows=1)
        assert dataset.shape[1] == 13
        variable = dataset[:, 1:].reshape(-1)
        timeseries = DataField(data=variable)
    else:
        raise ValueError('Unsupported dataset type: %s' % dataset_path)

    # cast datasat to DataField
    timeseries.create_time_array(
        date_from=datetime.strptime(FIRST_DATE, "%Y-%m-%d"), sampling='m')
    # assert dates as per npj paper
    timeseries.select_date(date(1900, 1, 1), date(2010, 1, 1))

    # prepare for surrogate creation
    if surrogates:
        seasonality = list(timeseries.get_seasonality(detrend=False))
        surrogate_field = SurrogateField()
        surrogate_field.copy_field(timeseries)

        timeseries.return_seasonality(seasonality[0], seasonality[1], None)
    print("Data loaded with shape %s" % (timeseries.data.shape))

    return timeseries, surrogate_field, seasonality


def compute_causality(timeseries1, timeseries2, tau_max, algorithm,
                      dim_condition=1, eta=0, phase_diff=False):
    """
    Compute causality as mean between 1 and tau_max.
    """
    results_simple = []
    results_knn = []
    for tau in range(1, tau_max):
        x, y, z = mutual.get_time_series_condition(
            [timeseries1, timeseries2], tau=tau,
            dim_of_condition=dim_condition, eta=eta, phase_diff=phase_diff)

        results_simple.append(mutual.cond_mutual_information(
            x, y, z, algorithm=algorithm, bins=NUM_BINS_EQQ))
        results_knn.append(mutual.knn_cond_mutual_information(
            x, y, z, k=K_KNN, dualtree=True))

    return np.mean(np.array(results_simple)), np.mean(np.array(results_knn))


def compute_information_measures(field, scales):
    """
    Computes information measures on a grid defined by scales on a given grid.
    :field: field to compute with
    :scales: scales for a grid
    """
    # prepare output arrays
    shape = (scales.shape[0], scales.shape[0])
    phase_phase_coherence = {'eqq': np.zeros(shape), 'knn': np.zeros(shape)}
    phase_amp_mi = {'eqq': np.zeros(shape), 'knn': np.zeros(shape)}
    phase_phase_causality = {'eqq': np.zeros(shape), 'knn': np.zeros(shape)}
    phase_amp_causality = {'eqq': np.zeros(shape), 'knn': np.zeros(shape)}

    for i, scale_i in enumerate(scales):
        field.wavelet(period=scale_i, period_unit='m', cut=1)
        phase_i = field.phase.copy()

        for j, scale_j in enumerate(scales):
            field.wavelet(period=scale_j, period_unit='m', cut=1)
            phase_j = field.phase.copy()
            amp_j = field.amplitude.copy()

            # compute phase-phase coherence
            phase_phase_coherence['eqq'][i, j] = mutual.mutual_information(
                phase_i, phase_j, algorithm='EQQ2', bins=NUM_BINS_EQQ)
            phase_phase_coherence['knn'][i, j] = mutual.knn_mutual_information(
                phase_i, phase_j, k=K_KNN, dualtree=True)

            # compute phase-amplitude mutual information
            phase_amp_mi['eqq'][i, j] = mutual.mutual_information(
                phase_i, amp_j, algorithm='EQQ2', bins=NUM_BINS_EQQ)
            phase_amp_mi['knn'][i, j] = mutual.knn_mutual_information(
                phase_i, amp_j, k=K_KNN, dualtree=True)

            # compute phase-phase causality
            eqq, knn = compute_causality(
                phase_i, phase_j, tau_max=7, algorithm="EQQ2", dim_condition=1,
                eta=0, phase_diff=True)
            phase_phase_causality['eqq'][i, j] = eqq
            phase_phase_causality['knn'][i, j] = knn

            # compute phase-amplitude causality
            eqq, knn = compute_causality(
                phase_i, np.power(amp_j, 2), tau_max=7, algorithm='GCM',
                dim_condition=3, eta=np.int(scale_i / 4), phase_diff=False)
            phase_amp_causality['eqq'][i, j] = eqq
            phase_amp_causality['knn'][i, j] = knn

    return (phase_phase_coherence, phase_amp_mi, phase_phase_causality,
            phase_amp_causality)


def _process_surrogates(field, seasonality, scales, jobq, resq):
    """
    Processes surrogates while job queue is not empty and puts results into
    result queue.
    """
    mean, var, _ = seasonality
    while True:
        s = jobq.get()
        # poison pill
        if s is None:
            break
        field.construct_fourier_surrogates(algorithm='FT')
        field.add_seasonality(mean, var, None)
        field.center_surr()
        surrogate_result = compute_information_measures(field, scales)
        resq.put(surrogate_result)


def main():
    # load timeseries
    timeseries, surrogate_field, seasonality = prepare_dataset(
        DATASET_PATH, "NINO3.4")
    timeseries.center_data()
    # get scales for computation
    scales = np.arange(PERIODS_SPAN[0], PERIODS_SPAN[-1] + 1, PERIODS_BIN)

    # compute for data
    print("Starting computing for data...")
    data_results = compute_information_measures(timeseries, scales)
    print("Data done!")
    data_results = ResultsContainer.from_algorithm_output(data_results)
    data_results.save(filename=SAVING_FILENAME + '_data.bin')

    # compute for surrogates
    print("Starting computing for %d surrogates using %d wokers" % (
        NUM_SURROGATES, WORKERS))
    surrogates_done = 0
    all_surrogates_results = []
    # prepare queues
    job_queue = Queue()
    result_queue = Queue()
    for _ in range(NUM_SURROGATES):
        job_queue.put(1)
    for _ in range(WORKERS):
        job_queue.put(None)

    workers = [Process(
        target=_process_surrogates, args=(surrogate_field, seasonality, scales,
                                          job_queue, result_queue))
               for i in range(WORKERS)]
    # start workers
    for worker in workers:
        worker.start()
    # fetch results
    progress_bar = tqdm(total=NUM_SURROGATES)
    while surrogates_done < NUM_SURROGATES:
        all_surrogates_results.append(result_queue.get())
        surrogates_done += 1
        progress_bar.update(1)
    for worker in workers:
        worker.join()
    print("Surrogates done.")
    surrogate_result = SurrogatesContainer.from_algorithm_output(
        all_surrogates_results)
    surrogate_result.save(filename=SAVING_FILENAME + '_surrogates.bin')


if __name__ == "__main__":
    main()
