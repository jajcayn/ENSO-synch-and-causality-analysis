import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Queue
import sys
sys.path.append('/home/nikola/Work/phd/multi-scale/src')
sys.path.append('/home/nikola/Work/phd/multi-scale/surrogates')

import wavelet_analysis
import mutual_information
from data_class import DataField
from surrogates import SurrogateField
import cPickle


def deseasonalize(ts):
    mean = np.zeros_like(ts)
    var = np.zeros_like(ts)
    mons = np.array([i%12 for i in range(ts.shape[0])])
    
    for m in range(12):
        sel = (m == mons)
        mean[sel] = np.mean(ts[sel], axis = 0)
        ts[sel] -= mean[sel]
        var[sel] = np.std(ts[sel], axis = 0, ddof = 1)
        ts[sel] /= var[sel]

    return mean, var


def phase_diff(ph1, ph2):
    ph = ph1 - ph2
    ph[ph < -np.pi] += 2*np.pi

    return ph


NEUTRAL = False
SPAN = [0.5, 7.5] # in years
NUM_SURR = 100
WRKRS = 4
ALG = 'EQQ2'

if NEUTRAL:
    enso = np.loadtxt("ensoPROneutral.txt")
else:
    enso = np.loadtxt("ensoPROdamped.txt")

# unit is month!!
enso = enso[:16384]

enso = DataField(data = enso)


enso_raw = enso.copy_data()
a = deseasonalize(enso.data)
enso_sg = SurrogateField(data = enso.data)


## annual phase
k0 = 6. # wavenumber of Morlet wavelet used in analysis
y = 12 # year in months
fourier_factor = (4 * np.pi) / (k0 + np.sqrt(2 + np.power(k0,2)))
period = y # frequency of interest
s0 = period / fourier_factor # get scale
wave, _, _, _ = wavelet_analysis.continous_wavelet(enso_raw, 1, False, wavelet_analysis.morlet, dj = 0, s0 = s0, j1 = 0, k0 = k0)
phaseAnn = np.arctan2(np.imag(wave), np.real(wave))[0, 1024:-1024]

## other scales
scales = np.arange(SPAN[0]*y, SPAN[1]*y+1, 1)
MI = []
CMI1 = []
CMI2 = []

print("Estimating MI // CMI on data...")

for sc in scales:
    s_temp = sc / fourier_factor # get scale
    wave, _, _, _ = wavelet_analysis.continous_wavelet(enso_raw, 1, False, wavelet_analysis.morlet, dj = 0, s0 = s_temp, j1 = 0, k0 = k0)
    phase_temp = np.arctan2(np.imag(wave), np.real(wave))[0, 1024:-1024]
 
    # MI
    MI.append(mutual_information.mutual_information(phase_temp, phaseAnn, algorithm = ALG, bins = 8))

    for tau in range(1, 6):
        condMI1 = []
        condMI2 = []

        condMI1.append(mutual_information.cond_mutual_information(phase_temp[:-tau], phase_diff(phaseAnn[tau:], phaseAnn[:-tau]), 
            phaseAnn[:-tau], algorithm = ALG, bins = 8))
        condMI2.append(mutual_information.cond_mutual_information(phaseAnn[:-tau], phase_diff(phase_temp[tau:], phase_temp[:-tau]), 
            phase_temp[:-tau], algorithm = ALG, bins = 8))

    CMI1.append(np.mean(np.array(condMI1)))
    CMI2.append(np.mean(np.array(condMI2)))

MI = np.array(MI)
CMI1 = np.array(CMI1)
CMI2 = np.array(CMI2)

print("Data done. Estimating on %d surrogates using %d workers..." % (NUM_SURR, WRKRS))


def _mi_surrs(sg, a, scales, phaseAnn, jobq, resq):
    mean, var = a
    while jobq.get() is not None:
        sg.construct_fourier_surrogates_spatial()
        sg.add_seasonality(mean, var, None)


        surrMI = []
        surrCMI1 = []
        surrCMI2 = []

        for sc in scales:
            s_temp = sc / fourier_factor # get scale
            wave, _, _, _ = wavelet_analysis.continous_wavelet(sg.surr_data, 1, False, wavelet_analysis.morlet, dj = 0, s0 = s_temp, j1 = 0, k0 = k0)
            phase_tempSurr = np.arctan2(np.imag(wave), np.real(wave))[0, 1024:-1024]

            surrMI.append(mutual_information.mutual_information(phase_tempSurr, phaseAnn, algorithm = 'EQQ', bins = 8))

            for tau in range(1, 6):
                condMI1 = []
                condMI2 = []

                condMI1.append(mutual_information.cond_mutual_information(phase_tempSurr[:-tau], phase_diff(phaseAnn[tau:], phaseAnn[:-tau]), 
                    phaseAnn[:-tau], algorithm = ALG, bins = 8))
                condMI2.append(mutual_information.cond_mutual_information(phaseAnn[:-tau], phase_diff(phase_tempSurr[tau:], phase_tempSurr[:-tau]), 
                    phase_tempSurr[:-tau], algorithm = ALG, bins = 8))

            surrCMI1.append(np.mean(np.array(condMI1)))
            surrCMI2.append(np.mean(np.array(condMI2)))

        resq.put((np.array(surrMI), np.array(surrCMI1), np.array(surrCMI2)))


surr_completed = 0
surrMI = np.zeros((NUM_SURR, MI.shape[0]))
surrCMI1 = np.zeros_like(surrMI)
surrCMI2 = np.zeros_like(surrMI)
jobq = Queue()
resq = Queue()
for i in range(NUM_SURR):
    jobq.put(1)
for i in range(WRKRS):
    jobq.put(None)

wrkrs = [Process(target=_mi_surrs, args = (enso_sg, a, scales, phaseAnn, jobq, resq))]
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

    fname = "zCMI_PRO%s_%dFTsurrs.bin" % ('neutral' if NEUTRAL else 'damped', NUM_SURR) 
    with open(fname, 'wb') as f:
        cPickle.dump({'zMI' : zMI, 'zCMI1' : zCMI1, 'zCMI2' : zCMI2}, f, protocol = cPickle.HIGHEST_PROTOCOL)

p = plt.subplot(311)
p.tick_params(axis='both', which='major', labelsize = 17)
p1, = plt.plot(scales, zMI, color = "#0059C7", linewidth = 2)
plt.xlim([scales[0], scales[-1]])
plt.ylim([-2, 10])
if NUM_SURR > 0:
 plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
plt.legend([p1], ['mutual information -- phaseAnn vs. phaseOthers'])
plt.xticks(scales[6::12], scales[6::12] / 12)
# plt.title("mutual information -- phaseAnn vs. phaseOthers")
p = plt.subplot(312)
p.tick_params(axis='both', which='major', labelsize = 17)
p1, = plt.plot(scales, zCMI1, color = "#0059C7", linewidth = 2)
plt.xlim([scales[0], scales[-1]])
plt.ylim([-2, 10])
plt.legend([p1], ['conditional mutual information -- phaseOthers -> phaseAnn'])
plt.xticks(scales[6::12], scales[6::12] / 12)
if NUM_SURR > 0:
    plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
    plt.ylabel("z-score", size = 20)
# plt.title("conditional mutual information -- phaseOthers -> phaseAnn")
p = plt.subplot(313)
p.tick_params(axis='both', which='major', labelsize = 17)
p1, = plt.plot(scales, zCMI2, color = "#0059C7", linewidth = 2)
plt.xlim([scales[0], scales[-1]])
plt.ylim([-2, 10])
if NUM_SURR > 0:
    plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
plt.legend([p1], ['conditional mutual information -- phaseAnn -> phaseOthers'])
plt.xticks(scales[6::12], scales[6::12] / 12)
plt.xlabel("period [years]", size = 20)
# plt.title("conditional mutual information -- phaseAnn -> phaseOthers")
plt.suptitle("PRO %s // z-score against %d FT surrogates" % ('neutral' if NEUTRAL else 'damped', NUM_SURR), size = 25)
plt.show()
# plt.savefig('test.png')