import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Queue
import sys
sys.path.append('/home/nikola/Work/phd/multi-scale/src')
sys.path.append('/home/nikola/Work/phd/multi-scale/surrogates')

import wavelet_analysis
import mutual_information as MI
from data_class import DataField
from surrogates import SurrogateField
import cPickle
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
from matplotlib.ticker import MultipleLocator, FuncFormatter

def load_enso_SSTs():
    # load enso SSTs
    print("[%s] Loading monthly ENSO SSTs..." % str(datetime.now()))
    enso_raw = np.loadtxt("nino34m13.txt") # length x 2 as 1st column is continuous year, second is SST in degC
    enso = DataField()

    enso.data = enso_raw[:, 1]

    time = np.zeros_like(enso.data, dtype = np.int32)
    y = np.int(enso_raw[0, 0])
    start_date = date(y, 1, 1)
    delta = relativedelta(months = +1)
    d = start_date
    for i in range(time.shape[0]):
        time[i] = d.toordinal()
        d += delta

    enso.time = time
    enso.location = 'NINO3.4 SSTs'

    # select 1024 data points
    enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = True)
    print("[%s] Data loaded with shape %s" % (str(datetime.now()), enso.data.shape))

    return enso

def phase_diff(ph1, ph2):
    ph = ph1 - ph2
    ph[ph < -np.pi] += 2*np.pi

    return ph

WVLT_SPAN = [5,93] # unit is month
NUM_SURR = 1000
WRKRS = 20
BINS = 4
COMPUTE = True # if True, the map will be evaluated, if False, it will be drawn

if COMPUTE:
    enso = load_enso_SSTs()

    ## DATA
    # prepare result matrices
    k0 = 6. # wavenumber of Morlet wavelet used in analysis
    y = 12 # year in months
    fourier_factor = (4 * np.pi) / (k0 + np.sqrt(2 + np.power(k0,2)))
    scales = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)
    phase_phase_coherence = np.zeros((scales.shape[0], scales.shape[0]))
    phase_phase_CMI = np.zeros_like(phase_phase_coherence)

    for i in range(phase_phase_coherence.shape[0]):
        sc_i = scales[i] / fourier_factor
        wave, _, _, _ = wavelet_analysis.continous_wavelet(enso.data, 1, False, wavelet_analysis.morlet, dj = 0, s0 = sc_i, j1 = 0, k0 = k0)
        phase_i = np.arctan2(np.imag(wave), np.real(wave))[0, 12:-12]
        
        for j in range(phase_phase_coherence.shape[1]):
            sc_j = scales[j] / fourier_factor
            wave, _, _, _ = wavelet_analysis.continous_wavelet(enso.data, 1, False, wavelet_analysis.morlet, dj = 0, s0 = sc_j, j1 = 0, k0 = k0)
            phase_j = np.arctan2(np.imag(wave), np.real(wave))[0, 12:-12]

            phase_phase_coherence[i, j] = MI.mutual_information(phase_i, phase_j, algorithm = 'EQQ2', bins = BINS)

            # conditional mutual inf -- avg over lags 1 - 6 months
            CMI = []
            for tau in range(1, 6):
                CMI.append(MI.cond_mutual_information(phase_i[:-tau], phase_diff(phase_j[tau:], phase_j[:-tau]), 
                    phase_j[:-tau], algorithm = 'EQQ2', bins = BINS))
            phase_phase_CMI[i, j] = np.mean(np.array(CMI))

    print("[%s] Analysis on data done." % str(datetime.now()))


    def _coh_cmi_surrs(sg, a, sc, jobq, resq):
        mean, var, _ = a
        while jobq.get() is not None:
            sg.construct_fourier_surrogates_spatial()
            sg.add_seasonality(mean, var, None)

            coh = np.zeros((sc.shape[0], sc.shape[0]))
            cmi = np.zeros_like(phase_phase_coherence)

            for i in range(coh.shape[0]):
                sc_i = sc[i] / fourier_factor
                wave, _, _, _ = wavelet_analysis.continous_wavelet(sg.surr_data, 1, False, wavelet_analysis.morlet, dj = 0, s0 = sc_i, j1 = 0, k0 = k0)
                phase_i = np.arctan2(np.imag(wave), np.real(wave))[0, 12:-12]
                
                for j in range(coh.shape[1]):
                    sc_j = sc[j] / fourier_factor
                    wave, _, _, _ = wavelet_analysis.continous_wavelet(sg.surr_data, 1, False, wavelet_analysis.morlet, dj = 0, s0 = sc_j, j1 = 0, k0 = k0)
                    phase_j = np.arctan2(np.imag(wave), np.real(wave))[0, 12:-12]

                    coh[i, j] = MI.mutual_information(phase_i, phase_j, algorithm = 'EQQ2', bins = BINS)

                    # conditional mutual inf -- avg over lags 1 - 6 months
                    CMI_temp = []
                    for tau in range(1, 6):
                        CMI_temp.append(MI.cond_mutual_information(phase_i[:-tau], phase_diff(phase_j[tau:], phase_j[:-tau]), 
                            phase_j[:-tau], algorithm = 'EQQ2', bins = BINS))
                    cmi[i, j] = np.mean(np.array(CMI_temp))

            resq.put((coh, cmi))


    ## SURROGATES
    if NUM_SURR > 0:
        print("[%s] Analysing %d FT surrogates using %d workers..." % (str(datetime.now()), NUM_SURR, WRKRS))
        a = enso.get_seasonality(DETREND = False)
        enso_sg = SurrogateField()
        enso_sg.copy_field(enso)

        surr_completed = 0
        surrCoherence = np.zeros(([NUM_SURR] + list(phase_phase_coherence.shape)))
        surrCMI = np.zeros_like(surrCoherence)
        jobq = Queue()
        resq = Queue()
        for i in range(NUM_SURR):
            jobq.put(1)
        for i in range(WRKRS):
            jobq.put(None)

        wrkrs = [Process(target=_coh_cmi_surrs, args = (enso_sg, a, scales, jobq, resq)) for i in range(WRKRS)]
        for w in wrkrs:
            w.start()

        while surr_completed < NUM_SURR:
            coh, cmi = resq.get()
            surrCoherence[surr_completed, :, :] = coh
            surrCMI[surr_completed, :, :] = cmi
            surr_completed += 1

            if surr_completed % 20 == 0:
                print("..%d/%d surrogate done.." % (surr_completed, NUM_SURR))

        for w in wrkrs:
            w.join()

        print("[%s] %d surrogates done. Saving..." % (str(datetime.now()), NUM_SURR))


    fname = ("CMImap%dbins.bin" % BINS)
    with open(fname, 'wb') as f:
        cPickle.dump({'phase x phase data' : phase_phase_coherence, 'phase CMI data' : phase_phase_CMI, 
            'phase x phase surrs' : surrCoherence, 'phase CMI surrs' : surrCMI}, f, protocol = cPickle.HIGHEST_PROTOCOL)

else:
    fname = ("CMImap%dbins.bin" % BINS)
    with open(fname, 'rb') as f:
        result = cPickle.load(f)

    phase_phase_coherence = result['phase x phase data']
    phase_phase_CMI = result['phase CMI data']
    surrCoherence = result['phase x phase surrs']
    surrCMI = result['phase CMI surrs']

    res_phase_coh = np.zeros_like(phase_phase_coherence)
    res_phase_cmi = np.zeros_like(res_phase_coh)

    for i in range(res_phase_coh.shape[0]):
        for j in range(res_phase_coh.shape[1]):
            res_phase_coh[i, j] = np.sum(np.greater(phase_phase_coherence[i, j], surrCoherence[:, i, j])) / np.float(surrCoherence.shape[0])
            res_phase_cmi[i, j] = np.sum(np.greater(phase_phase_CMI[i, j], surrCMI[:, i, j])) / np.float(surrCMI.shape[0])

    scales = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)
    # a = np.random.rand(scales.shape[0], scales.shape[0]) + 0.5
    x, y = np.meshgrid(scales, scales)
    fig, axs = plt.subplots(1, 2, figsize = (13,7))
    i = 1
    for ax in axs.ravel():
        cs = ax.contourf(x, y, a, levels = np.arange(0.95, 1, 0.00125), cmap = plt.cm.get_cmap("jet"))
        ax.tick_params(axis='both', which='major', labelsize = 17)
        ax.set_title("test")
        ax.xaxis.set_major_locator(MultipleLocator(12))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
        ax.xaxis.set_minor_locator(MultipleLocator(6))
        ax.yaxis.set_major_locator(MultipleLocator(12))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
        ax.yaxis.set_minor_locator(MultipleLocator(6))
        ax.set_xlabel("period [years]", size = 20)
        if i == 1:
            ax.set_ylabel("period [years]", size = 20)
        i += 1
        # fig.colorbar(cs, ax = ax, shrink = 0.8)

    # plt.savefig('test.png')
    print scales



