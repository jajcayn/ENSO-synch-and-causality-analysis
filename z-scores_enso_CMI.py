import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Queue
import sys
sys.path.append('/home/nikola/Work/multi-scale')
# sys.path.append('/Users/nikola/work-ui/multi-scale')
import src.wavelet_analysis as wvlt
import src.mutual_information as mutual_information
from src.data_class import DataField
import cPickle
from datetime import datetime, date
from dateutil.relativedelta import relativedelta


def load_enso_SSTs(emr = False, num_ts = 0):
    from src.data_class import DataField
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

    print enso.data.shape
    print enso.get_date_from_ndx(0), enso.get_date_from_ndx(-1)

    # select 1024 data points
    enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = True)

    if emr:
        import scipy.io as sio
        print("using synthetic time series from ERM number %d" % (num_ts+1))
        raw = sio.loadmat("Nino34-ERM-1884-2013linear-same-init-cond-no-anomalies.mat")['N34s']
        enso.data = raw[-1024:, num_ts]

    print("[%s] Data loaded with shape %s" % (str(datetime.now()), enso.data.shape))

    return enso


def phase_diff(ph1, ph2):
    ph = ph1 - ph2
    ph[ph < -np.pi] += 2*np.pi

    return ph


def get_seasonality_arbitrary(ts):
    """
    assumes monthly data
    outputs as DataField.get_seasonality() function
    """

    mean = np.zeros_like(ts)
    std = np.zeros_like(ts)

    for i in range(12):
        sel = slice(i, ts.shape[0], 12)
        mean[sel] = np.nanmean(ts[sel], axis = 0)
        ts[sel] -= mean[sel]

        std[sel] = np.nanstd(ts[sel], axis = 0, ddof = 1)
        std[std == 0.] = 1.
        ts[sel] /= std[sel]

    return ts, (mean, std, None)


def _mi_surrs(sg, a, scales, phaseAnn, jobq, resq):
    mean, var, _ = a
    while jobq.get() is not None:
        if "nino34" in DATA:
            sg.construct_fourier_surrogates_spatial()
            sg.add_seasonality(mean, var, None)
            surrogate = sg.surr_data.copy()
        elif DATA == "PRO":
            surrogate = get_single_FT_surrogate(sg)
            surrogate *= var
            surrogate += mean

        surrMI = []
        surrCMI1 = []
        surrCMI2 = []

        for sc in scales:
            s_temp = sc / fourier_factor # get scale
            wave, _, _, _ = wvlt.continous_wavelet(surrogate, 1, False, wvlt.morlet, dj = 0, s0 = s_temp, j1 = 0, k0 = k0)
            phase_tempSurr = np.arctan2(np.imag(wave), np.real(wave))[0, 12:-12]

            surrMI.append(mutual_information.mutual_information(phase_tempSurr, phaseAnn, algorithm = 'EQQ', bins = BINS))

            for tau in range(1, 6):
                condMI1 = []
                condMI2 = []

                condMI1.append(mutual_information.cond_mutual_information(phase_tempSurr[:-tau], phase_diff(phaseAnn[tau:], phaseAnn[:-tau]), 
                    phaseAnn[:-tau], algorithm = ALG, bins = BINS))
                condMI2.append(mutual_information.cond_mutual_information(phaseAnn[:-tau], phase_diff(phase_tempSurr[tau:], phase_tempSurr[:-tau]), 
                    phase_tempSurr[:-tau], algorithm = ALG, bins = BINS))

            surrCMI1.append(np.mean(np.array(condMI1)))
            surrCMI2.append(np.mean(np.array(condMI2)))

        resq.put((np.array(surrMI), np.array(surrCMI1), np.array(surrCMI2)))



DATA = "nino34ERM" # "nino34" or "PRO"
SPAN = [0.5, 7.5] # in years
NUM_SURR = 100
WRKRS = 7
ALG = 'EQQ2'
BINS = 4

for model_no in range(100):

    if DATA == "nino34":
        enso = load_enso_SSTs()
    elif DATA == "nino34ERM":
        enso = load_enso_SSTs(emr = True, num_ts = model_no)
    elif DATA == "PRO":
        from parametric_recharge_oscillator import ENSO_PROmodel
        enso = ENSO_PROmodel(length = 1024, daily = False, damped = True, ENSOperiod = 3.75, modulation = 2, lambda0 = 0.4)
        enso.integrate_PROmodel()


    ## annual phase
    k0 = 6. # wavenumber of Morlet wavelet used in analysis
    y = 12 # year in months
    fourier_factor = (4 * np.pi) / (k0 + np.sqrt(2 + np.power(k0,2)))
    period = y # frequency of interest
    s0 = period / fourier_factor # get scale
    wave, _, _, _ = wvlt.continous_wavelet(enso.data, 1, False, wvlt.morlet, dj = 0, s0 = s0, j1 = 0, k0 = k0)
    phaseAnn = np.arctan2(np.imag(wave), np.real(wave))[0, 12:-12]

    ## other scales
    scales = np.arange(SPAN[0]*y, SPAN[1]*y+1, 1)
    MI = []
    CMI1 = []
    CMI2 = []

    print("Estimating MI // CMI on data...")

    for sc in scales:
        s_temp = sc / fourier_factor # get scale
        wave, _, _, _ = wvlt.continous_wavelet(enso.data, 1, False, wvlt.morlet, dj = 0, s0 = s_temp, j1 = 0, k0 = k0)
        phase_temp = np.arctan2(np.imag(wave), np.real(wave))[0, 12:-12]
     
        # MI
        MI.append(mutual_information.mutual_information(phase_temp, phaseAnn, algorithm = ALG, bins = BINS))

        for tau in range(1, 6):
            condMI1 = []
            condMI2 = []

            condMI1.append(mutual_information.cond_mutual_information(phase_temp[:-tau], phase_diff(phaseAnn[tau:], phaseAnn[:-tau]), 
                phaseAnn[:-tau], algorithm = ALG, bins = BINS))
            condMI2.append(mutual_information.cond_mutual_information(phaseAnn[:-tau], phase_diff(phase_temp[tau:], phase_temp[:-tau]), 
                phase_temp[:-tau], algorithm = ALG, bins = BINS))

        CMI1.append(np.mean(np.array(condMI1)))
        CMI2.append(np.mean(np.array(condMI2)))

    MI = np.array(MI)
    CMI1 = np.array(CMI1)
    CMI2 = np.array(CMI2)


    print("Data done. Estimating on %d surrogates using %d workers..." % (NUM_SURR, WRKRS))

    if NUM_SURR > 0:
        surr_completed = 0
        from src.surrogates import SurrogateField, get_single_FT_surrogate

        if DATA == "PRO":
            enso_sg, a = get_seasonality_arbitrary(enso.data)
        elif "nino34" in DATA:
            a = enso.get_seasonality(DETREND = False)
            enso_sg = SurrogateField()
            enso_sg.copy_field(enso)


        surrMI = np.zeros((NUM_SURR, MI.shape[0]))
        surrCMI1 = np.zeros_like(surrMI)
        surrCMI2 = np.zeros_like(surrMI)
        jobq = Queue()
        resq = Queue()
        for i in range(NUM_SURR):
            jobq.put(1)
        for i in range(WRKRS):
            jobq.put(None)

        wrkrs = [Process(target=_mi_surrs, args = (enso_sg, a, scales, phaseAnn, jobq, resq)) for i in range(WRKRS)]
        for w in wrkrs:
            w.start()

        while surr_completed < NUM_SURR:
            mi, cmi1, cmi2 = resq.get()
            surrMI[surr_completed, :] = mi
            surrCMI1[surr_completed, :] = cmi1
            surrCMI2[surr_completed, :] = cmi2
            surr_completed += 1

            # if surr_completed % 20 == 0:
            print("..%d/%d surrogate done.." % (surr_completed, NUM_SURR))

        for w in wrkrs:
            w.join()

        print("Surrogates done. Plotting the z-score..")

        MIstd = np.std(surrMI, axis = 0, ddof = 1)
        MIstd[MIstd == 0] = 1.
        CMI1std = np.std(surrCMI1, axis = 0, ddof = 1)
        CMI1std[CMI1std == 0] = 1.
        CMI2std = np.std(surrCMI2, axis = 0, ddof = 1)
        CMI2std[CMI2std == 0] = 1.

        if NUM_SURR > 0:
            zMI = (MI - np.mean(surrMI, axis = 0)) / MIstd
            zCMI1 = (CMI1 - np.mean(surrCMI1, axis = 0)) / CMI1std
            zCMI2 = (CMI2 - np.mean(surrCMI2, axis = 0)) / CMI2std

            fname = "bins/zCMI_%sSST_%dFTsurrs_-SIC-no-anom-TSno%d.bin" % (DATA.upper(), NUM_SURR, model_no) 
            with open(fname, 'wb') as f:
                cPickle.dump({'zMI' : zMI, 'zCMI1' : zCMI1, 'zCMI2' : zCMI2,
                    'MI' : MI, 'CMI1' : CMI1, 'CMI2' : CMI2,
                    'surrMI' : surrMI, 'surrCMI1' : surrCMI1, 'surrCMI2' : surrCMI2}, 
                    f, protocol = cPickle.HIGHEST_PROTOCOL)

        plt.figure()
        p = plt.subplot(311)
        p.tick_params(axis='both', which='major', labelsize = 17)
        p1, = plt.plot(scales, zMI, color = "#0059C7", linewidth = 2)
        # print zMI.shape
        plt.xlim([scales[0], scales[-1]])
        plt.ylim([-2, 11])
        if NUM_SURR > 0:
            plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
        plt.legend([p1], ['mutual information -- phaseAnn vs. phaseOthers'])
        plt.xticks(scales[6::12], scales[6::12] / 12)
        # plt.title("mutual information -- phaseAnn vs. phaseOthers")
        p = plt.subplot(313)
        p.tick_params(axis='both', which='major', labelsize = 17)
        p1, = plt.plot(scales, zCMI1, color = "#0059C7", linewidth = 2)
        plt.xlim([scales[0], scales[-1]])
        plt.ylim([-2, 18])
        plt.legend([p1], ['conditional mutual information -- phaseOthers -> phaseAnn'])
        plt.xticks(scales[6::12], scales[6::12] / 12)
        if NUM_SURR > 0:
            plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
        plt.xlabel("period [years]", size = 20)
        # plt.title("conditional mutual information -- phaseOthers -> phaseAnn")
        p = plt.subplot(312)
        p.tick_params(axis='both', which='major', labelsize = 17)
        p1, = plt.plot(scales, zCMI2, color = "#0059C7", linewidth = 2)
        plt.xlim([scales[0], scales[-1]])
        plt.ylim([-2, 8])
        plt.legend([p1], ['conditional mutual information -- phaseAnn -> phaseOthers'])
        plt.xticks(scales[6::12], scales[6::12] / 12)
        if NUM_SURR > 0:
            plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
            plt.ylabel("z-score", size = 20)
        # plt.title("conditional mutual information -- phaseAnn -> phaseOthers")
        if DATA == "nino34":
            plt.suptitle("NINO3.4 SST // z-score against %d FT surrogates" % (NUM_SURR), size = 25)
        elif DATA == "nino34ERM":
            plt.suptitle("ERM model 1884-2013 linear %d. ts // z-score against %d FT surrogates" % (model_no+1, NUM_SURR), size = 25)
        elif DATA == "PRO":
            plt.suptitle("PRO model -- damped SST // z-score against %d FT surrogates" % (NUM_SURR), size = 25)
        # plt.show()
        plt.savefig('plots/ERM1884-2013linear-SIC-no-anom_z-score_annual%d.png' % (model_no))
        plt.close()