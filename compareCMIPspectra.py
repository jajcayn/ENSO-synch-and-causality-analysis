import numpy as np
import cPickle


def get_L_distance(a, b, L = 1):
    res = 0
    if a.shape[0] != b.shape[0]:
        raise Exception("time series must have same length")
    for t in range(a.shape[0]):
        res += np.power(np.abs(a[t] - b[t]), L)
    return np.power(res, 1./L)

CMIP5models = ['CanESM2', 'GFDLCM3', 'GISSE2Hp1', 'GISSE2Hp2', 'GISSE2Hp3', 'GISSE2Rp1']
CMIP5models += ['GISSE2Rp2', 'GISSE2Rp3', 'HadGem2ES', 'IPSL_CM5A_LR', 'MIROC5', 'MRICGCM3']
CMIP5models += ['CCSM4', 'CNRMCM5', 'CSIROmk360']


with open("CMIP5spectra/nino34.bin", "rb") as f:
    data = cPickle.load(f)
    scales = data['scales']
    nino_pow = data['wvlt_power']
    nino_coh = data['autocoherence_ph']

results = []
models = []
for CMIP5model in CMIP5models:
    for model_count in range(10):
        try:
            with open("CMIP5spectra/SPECTRUM-%s-%d.bin" % (CMIP5model, model_count), "rb") as f:
                data = cPickle.load(f)
            power = data['wvlt_power']
            coh = data['autocoherence_ph']
        except IOError:
            continue

        temp = {}
        temp['wvltL1'] = get_L_distance(nino_pow, power, 1)
        temp['wvltL2'] = get_L_distance(nino_pow, power, 2)
        temp['cohL1'] = get_L_distance(nino_coh, coh, 1)
        temp['cohL2'] = get_L_distance(nino_coh, coh, 2)

        models.append(CMIP5model + "_" + str(model_count))
        results.append(temp)


# items = sorted(results, key=lambda x: x['wvltL1'])
ndxwL1 = [i[0] for i in sorted(enumerate(results), key=lambda x:x[1]['wvltL1'])]
ndxwL2 = [i[0] for i in sorted(enumerate(results), key=lambda x:x[1]['wvltL2'])]
ndxcL1 = [i[0] for i in sorted(enumerate(results), key=lambda x:x[1]['cohL1'])]
ndxcL2 = [i[0] for i in sorted(enumerate(results), key=lambda x:x[1]['cohL2'])]

import csv
with open("spectraCMIPresults.csv", "w") as f:
    writer = csv.writer(f)

    writer.writerow(["BEST WVLT L1", "BEST WVLT L1 value", "BEST WVLT L2", "BEST WVLT L2 value", "BEST COH L1", "BEST COH L1 value", "BEST COH L2", "BEST COH L2 value"])
    for a, b, c, d in zip(ndxwL1, ndxwL2, ndxcL1, ndxcL2):
        writer.writerow([models[a], results[a]['wvltL1'], models[b], results[b]['wvltL2'], models[c], results[c]['cohL1'], models[d], results[d]['cohL2']])
# for l in ndxwL1:
#     print models[l], results[l]['wvltL1']





