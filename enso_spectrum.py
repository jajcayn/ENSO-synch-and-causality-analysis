import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
from matplotlib.ticker import MultipleLocator, FuncFormatter
import cPickle
import sys
import matplotlib.gridspec as gridspec
sys.path.append('/Users/nikola/work-ui/multi-scale')
import src.wavelet_analysis as wvlt
import src.mutual_information as MI
from src.data_class import DataField
from src.surrogates import SurrogateField



def load_enso_SSTs(CMIP5model = None, num_ts = None, PROmodel = False):
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

    if CMIP5model is not None and num_ts is not None:
        print("Loading %s model, time series %d..." % (CMIP5model, num_ts))
        fname = CMIP5model + '.txt'
        model = np.loadtxt('N34_CMIP5/' + fname)
        enso.data = model[:, num_ts]

        # varnorm, not seasonal
        enso.data /= np.std(enso.data, axis = 0, ddof = 1)

    if PROmodel:
        print("[%s] Integrating PRO model which will be used instead of ENSO SSTs..." % (str(datetime.now())))
        from parametric_recharge_oscillator import ENSO_PROmodel
        PROmodel_enso = ENSO_PROmodel(length = enso.data.shape[0], daily = False, damped = True, ENSOperiod = 3.75, modulation = 2, lambda0 = 0.4)
        PROmodel_enso.integrate_PROmodel()
        enso.data = PROmodel_enso.data.copy()

    if NUM_SURR > 0:
        a = list(enso.get_seasonality(DETREND = False))
        enso_sg = SurrogateField()

        _, _, idx = enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = False)
        enso_sg.copy_field(enso)
        enso_sg.data = enso_sg.data[idx[0]:idx[1]]

        enso.return_seasonality(a[0], a[1], None)

        a[0] = a[0][idx[0]:idx[1]]
        a[1] = a[1][idx[0]:idx[1]]

    # select 1024 data points
    enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = True)
    print("[%s] Data loaded with shape %s" % (str(datetime.now()), enso.data.shape))

    return enso, enso_sg, a


# CMIP5models = ['N34_CanESM2', 'N34_GFDLCM3', 'N34_GISSE2Hp1', 'N34_GISSE2Hp2', 'N34_GISSE2Hp3', 'N34_GISSE2Rp1']
# CMIP5models += ['N34_GISSE2Rp2', 'N34_GISSE2Rp3', 'N34_HadGem2ES', 'N34_IPSL_CM5A_LR', 'N34_MIROC5', 'N34_MRICGCM3']
# CMIP5models += ['N34_CCSM4', 'N34_CNRMCM5', 'N34_CSIROmk360']

NUM_SURR = 1000
use_PRO_model = False
CMIP5models = [None]


for CMIP5model in CMIP5models:
    # fname = CMIP5model + '.txt'
    # model = np.loadtxt('N34_CMIP5/' + fname)
    # model_count = model.shape[1]
    import scipy.io as sio
    a = sio.loadmat("Nino34-ERM-1884-2013linear-same-init-cond-no-anomalies.mat")
    sim_nino = a['N34s'] # 1920 x 100 as ts length x ensemble
    model_count = sim_nino.shape[1]
    sim_nino = sim_nino[-1024:, :]

    plt1 = []
    plt2 = []
    plt3 = []
    plt4 = []
    plt5 = []
    plt6 = []

    # data_type = "PRO model damped integrated monthly data -- 100members"
    # fname = "SPECTRUM-PROdamped-ensemble100.png"

    for num_ts in range(model_count):
        print num_ts

        # data_type = ("CMIP5 -- %s / %d monthly data" % (CMIP5model[4:], num_ts))
        # fname = ("SPECTRUM-%s-%d.bin" % (CMIP5model[4:], num_ts))

        # enso, enso_sg, seasonality = load_enso_SSTs(CMIP5model = CMIP5model, num_ts = num_ts, PROmodel = False)
        # enso.data = sim_nino[:, num_ts].copy()
        k0 = 6. # wavenumber of Morlet wavelet used in analysis
        y = 12 # year in months
        fourier_factor = (4 * np.pi) / (k0 + np.sqrt(2 + np.power(k0,2)))
        WVLT_SPAN = [5,96] # in months
        scales = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)

        wvlt_power = np.zeros((scales.shape[0],))
        autocoherence_ph = np.zeros((scales.shape[0],))
        autocoherence_re = np.zeros_like(autocoherence_ph)

        for sc, i in zip(scales, range(wvlt_power.shape[0])):
            # print sc
            scale = sc / fourier_factor
            wave, _, _, _ = wvlt.continous_wavelet(sim_nino[:, num_ts], 1, False, wvlt.morlet, dj = 0, s0 = scale, j1 = 0, k0 = k0)
            wave = wave[0, 12:-12]
            phase = np.arctan2(np.imag(wave), np.real(wave))

            # wvlt power
            wvlt_power[i] = np.sum(np.power(np.abs(wave), 2))

            # autocoherence
            autocoh = []
            autocoh2 = []
            for tau in range(1, 90):
                x, y, _ = MI.get_time_series_condition([phase, phase], tau = tau, dim_of_condition = 0)
                autocoh.append(MI.mutual_information(x, y, algorithm = 'EQQ2', bins = 4, log2 = False))
                x, y, _ = MI.get_time_series_condition([np.real(wave), np.real(wave)], tau = tau, dim_of_condition = 0)
                autocoh2.append(MI.mutual_information(x, y, algorithm = 'EQQ2', bins = 4, log2 = False))
            autocoherence_ph[i] = np.mean(autocoh)
            autocoherence_re[i] = np.mean(autocoh2)


        WVLT_SPAN = [90, 360]
        scales2 = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)

        wvlt_power2 = np.zeros((scales2.shape[0],))
        autocoherence_ph2 = np.zeros((scales2.shape[0],))
        autocoherence_re2 = np.zeros_like(autocoherence_ph2)

        for sc, i in zip(scales2, range(wvlt_power2.shape[0])):
            # print sc
            scale = sc / fourier_factor
            wave, _, _, _ = wvlt.continous_wavelet(sim_nino[:, num_ts], 1, False, wvlt.morlet, dj = 0, s0 = scale, j1 = 0, k0 = k0)
            wave = wave[0, 12:-12]
            phase = np.arctan2(np.imag(wave), np.real(wave))

            # wvlt power
            wvlt_power2[i] = np.sum(np.power(np.abs(wave), 2))

            # autocoherence
            autocoh = []
            autocoh2 = []
            for tau in range(1, 90):
                x, y, _ = MI.get_time_series_condition([phase, phase], tau = tau, dim_of_condition = 0)
                autocoh.append(MI.mutual_information(x, y, algorithm = 'EQQ2', bins = 4, log2 = False))
                x, y, _ = MI.get_time_series_condition([np.real(wave), np.real(wave)], tau = tau, dim_of_condition = 0)
                autocoh2.append(MI.mutual_information(x, y, algorithm = 'EQQ2', bins = 4, log2 = False))
            autocoherence_ph2[i] = np.mean(autocoh)
            autocoherence_re2[i] = np.mean(autocoh2)

        wvlt_power /= 1000
        wvlt_power2 /= 1000
        plt1.append(autocoherence_ph)
        plt2.append(autocoherence_re)
        plt3.append(wvlt_power)
        plt4.append(autocoherence_ph2)
        plt5.append(autocoherence_re2)
        plt6.append(wvlt_power2)

        # with open("spectra/" + fname, "wb") as f:
        #     cPickle.dump({"scales" : scales, "autocoherence_re" : autocoherence_re, "autocoherence_ph" : autocoherence_ph, 
        #         "wvlt_power" : wvlt_power, "scales2" : scales2, "autocoherence_re2" : autocoherence_re2, "autocoherence_ph2" : autocoherence_ph2, 
        #         "wvlt_power2" : wvlt_power2}, f, protocol = cPickle.HIGHEST_PROTOCOL)


    plt1 = np.array(plt1)
    plt2 = np.array(plt2)
    plt3 = np.array(plt3)
    plt4 = np.array(plt4)
    plt5 = np.array(plt5)
    plt6 = np.array(plt6)


    fig = plt.figure(figsize=(20,15))
    gs = gridspec.GridSpec(2, 2)
    gs.update(left=0.08, right=0.92, hspace=0.05, top=0.93, bottom=0.05, wspace=0.12)
    ax = plt.subplot(gs[0,0])
    ax.set_xlim(scales[0], scales[-1])
    ax.xaxis.set_major_locator(MultipleLocator(12))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
    ax.xaxis.set_minor_locator(MultipleLocator(6))
    ax.tick_params(axis='both', which='major', labelsize = 22)
    ax.tick_params(axis='x', which='both', labelbottom='off')
    for vl in range(6,scales[-1],6):
        ax.axvline(vl, color = "#9B9B9B")

    # ax.plot(scales, autocoherence_ph, linewidth = 4, color = "#3C0B1F", label = "phase $\phi$")
    # ax.plot(scales, autocoherence_re, linewidth = 4, color = "#F67A48", label = "Re(W(t))")
    ax.errorbar(scales, np.mean(plt1, axis = 0), yerr = np.std(plt1, ddof = 1, axis = 0), linewidth = 4, color = "#3C0B1F", label = "phase $\phi$")
    ax.errorbar(scales, np.mean(plt2, axis = 0), yerr = np.std(plt2, ddof = 1, axis = 0), linewidth = 4, color = "#F67A48", label = "Re(W(t))")
    ax.legend()

    ax.set_ylabel("AUTOCOHERENCE", size = 27)

    ax = plt.subplot(gs[1,0])
    ax.set_xlim(scales[0], scales[-1])
    ax.xaxis.set_major_locator(MultipleLocator(12))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
    ax.xaxis.set_minor_locator(MultipleLocator(6))
    ax.tick_params(axis='both', which='major', labelsize = 22)
    for vl in range(6,scales[-1],6):
        ax.axvline(vl, color = "#9B9B9B")

    ax.set_ylabel("WAVELET POWER [*1000]", size = 27)
    ax.set_xlabel("PERIOD [years]", size = 27)
    # ax.plot(scales, wvlt_power, linewidth = 4, color = "#3C0B1F", label = "phase $\phi$")
    ax.errorbar(scales, np.mean(plt3, axis = 0), yerr = np.std(plt3, ddof = 1, axis = 0), linewidth = 4, color = "#3C0B1F")


    ax = plt.subplot(gs[0,1])
    ax.set_xlim(scales2[0], scales2[-1])
    ax.xaxis.set_major_locator(MultipleLocator(24))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
    ax.xaxis.set_minor_locator(MultipleLocator(6))
    ax.tick_params(axis='both', which='major', labelsize = 22)
    ax.tick_params(axis='x', which='both', labelbottom='off')
    for vl in range(12,scales2[-1],12):
        ax.axvline(vl, color = "#9B9B9B")
    # ax.plot(scales2, autocoherence_ph2, linewidth = 4, color = "#3C0B1F", label = "phase $\phi$")
    # ax.plot(scales2, autocoherence_re2, linewidth = 4, color = "#F67A48", label = "Re(W(t))")
    ax.errorbar(scales2, np.mean(plt4, axis = 0), yerr = np.std(plt4, ddof = 1, axis = 0), linewidth = 4, color = "#3C0B1F", label = "phase $\phi$")
    ax.errorbar(scales2, np.mean(plt5, axis = 0), yerr = np.std(plt5, ddof = 1, axis = 0), linewidth = 4, color = "#F67A48", label = "Re(W(t))")
    ax.legend()

    # ax.set_ylabel(r"AUTOCOHERENCE $I(x(t); x(t+\tau))$ [nats]", size = 27)

    ax = plt.subplot(gs[1,1])
    ax.set_xlim(scales2[0], scales2[-1])
    ax.xaxis.set_major_locator(MultipleLocator(24))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
    ax.xaxis.set_minor_locator(MultipleLocator(6))
    ax.tick_params(axis='both', which='major', labelsize = 22)
    for vl in range(12,scales2[-1],12):
        ax.axvline(vl, color = "#9B9B9B")

    # ax.set_ylabel("WAVELET POWER", size = 27)
    ax.set_xlabel("PERIOD [years]", size = 27)
    # ax.plot(scales2, wvlt_power2, linewidth = 4, color = "#3C0B1F", label = "phase $\phi$")
    ax.errorbar(scales2, np.mean(plt6, axis = 0), yerr = np.std(plt6, ddof = 1, axis = 0), linewidth = 4, color = "#3C0B1F")

    plt.suptitle("Regression model -- 100 realisations", size = 35)

    # plt.savefig("spectra/" + fname)
    plt.savefig("spectra/simulatedNINO34.png")

