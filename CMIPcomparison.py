import numpy as np
import cPickle
import csv
from sklearn.metrics import adjusted_rand_score

def get_L_distance(a, b, L = 1):
    res = 0
    if a.shape[0] != b.shape[0]:
        raise Exception("time series must have same length")
    for t in range(a.shape[0]):
        res += np.power(np.abs(a[t] - b[t]), L)
    return np.power(res, 1./L)


def evaluate_MI(fname, threshold = 0.95):
    CUT = slice(0,1000)
    # version = 3
    with open(fname, 'rb') as f:
        result = cPickle.load(f)

    phase_phase_coherence = result['phase x phase data']
    phase_phase_CMI = result['phase CMI data']
    surrCoherence = result['phase x phase surrs'][CUT, ...]
    surrCMI = result['phase CMI surrs'][CUT, ...]
    phase_amp_condMI = result['phase amp CMI data']
    surrPhaseAmpCMI = result['phase amp CMI surrs'][CUT, ...]

    res_phase_coh = np.zeros_like(phase_phase_coherence)
    res_phase_cmi = np.zeros_like(res_phase_coh)
    res_phase_amp_CMI = np.zeros_like(res_phase_coh)

    for i in range(res_phase_coh.shape[0]):
        for j in range(res_phase_coh.shape[1]):
            res_phase_coh[i, j] = np.sum(np.greater(phase_phase_coherence[i, j], surrCoherence[:, i, j])) / np.float(surrCoherence.shape[0])
            res_phase_cmi[i, j] = np.sum(np.greater(phase_phase_CMI[i, j], surrCMI[:, i, j])) / np.float(surrCMI.shape[0])
            res_phase_amp_CMI[i, j] = np.sum(np.greater(phase_amp_condMI[i, j], surrPhaseAmpCMI[:, i, j])) / np.float(surrPhaseAmpCMI.shape[0])

    f.close()

    res_phase_coh_thr = np.zeros_like(res_phase_coh, dtype = np.int)
    res_phase_coh_thr[np.where(res_phase_coh > threshold)] = 1
    res_phase_cmi_thr = np.zeros_like(res_phase_cmi, dtype = np.int)
    res_phase_cmi_thr[np.where(res_phase_cmi > threshold)] = 1
    res_phase_amp_CMI_thr = np.zeros_like(res_phase_amp_CMI, dtype = np.int)
    res_phase_amp_CMI_thr[np.where(res_phase_amp_CMI > threshold)] = 1

    return res_phase_coh_thr, res_phase_cmi_thr, res_phase_amp_CMI_thr

WVLT_SPAN = [5,93]
CMIP5models = ['CanESM2', 'GFDLCM3', 'GISSE2Hp1', 'GISSE2Hp2', 'GISSE2Hp3', 'GISSE2Rp1']
CMIP5models += ['GISSE2Rp2', 'GISSE2Rp3', 'HadGem2ES', 'IPSL_CM5A_LR', 'MIROC5', 'MRICGCM3']
CMIP5models += ['CCSM4', 'CNRMCM5', 'CSIROmk360']


with open("CMIP5spectra/nino34.bin", "rb") as f:
    data = cPickle.load(f)
    scales = data['scales']
    nino_pow = data['wvlt_power']

fname = ("results/KNN_CMImap_k_32_3Dcond_GaussCorr.bin")
nino34_inter = evaluate_MI(fname)

results = []
models = []
total = 0
modelscount = 0
for CMIP5model in CMIP5models:
    modelscount += 1
    for model_count in range(10):
        try:
            with open("CMIP5spectra/SPECTRUM-%s-%d.bin" % (CMIP5model, model_count), "rb") as f:
                data = cPickle.load(f)
            power = data['wvlt_power']
            coh = data['autocoherence_ph']
            result_temp = evaluate_MI("models/kNN_CMImap_k_32_3DcondN34_%sts%d.bin" % (CMIP5model, model_count))
            total += 1
        except IOError:
            break

        temp = {}
        # spectra
        temp['wvltL1'] = get_L_distance(nino_pow, power, 1)
        temp['wvltL2'] = get_L_distance(nino_pow, power, 2)
        temp['corr'] = np.corrcoef(nino_pow, power)[0, 1]

        # interactions
        for i, metric in zip(range(len(result_temp)), ['synch', 'ph-caus', 'ph-amp-caus']):
            temp['%s-corr' % metric] = np.corrcoef(nino34_inter[i].reshape(-1), result_temp[i].reshape(-1))[0,1]
            temp[metric+'-ARI'] = adjusted_rand_score(nino34_inter[i].reshape(-1), result_temp[i].reshape(-1))

        models.append(CMIP5model + "_" + str(model_count))
        results.append(temp)

    models.append(CMIP5model + "_ens_mean")
    last = results[-model_count:]
    temp = {}
    temp['wvltL1'] = np.mean([a['wvltL1'] for a in last])
    temp['wvltL2'] = np.mean([a['wvltL2'] for a in last])
    temp['corr'] = np.mean([a['corr'] for a in last])
    temp['synch-corr'] = np.mean([a['synch-corr'] for a in last])
    temp['synch-ARI'] = np.mean([a['synch-ARI'] for a in last])
    temp['ph-caus-corr'] = np.mean([a['ph-caus-corr'] for a in last])
    temp['ph-caus-ARI'] = np.mean([a['ph-caus-ARI'] for a in last])
    temp['ph-amp-caus-corr'] = np.mean([a['ph-amp-caus-corr'] for a in last])
    temp['ph-amp-caus-ARI'] = np.mean([a['ph-amp-caus-ARI'] for a in last])
    results.append(temp)

    models.append(CMIP5model + "_ens_std")
    temp = {}
    temp['wvltL1'] = np.std([a['wvltL1'] for a in last], ddof = 1)
    temp['wvltL2'] = np.std([a['wvltL2'] for a in last], ddof = 1)
    temp['corr'] = np.std([a['corr'] for a in last], ddof = 1)
    temp['synch-corr'] = np.std([a['synch-corr'] for a in last], ddof = 1)
    temp['synch-ARI'] = np.std([a['synch-ARI'] for a in last], ddof = 1)
    temp['ph-caus-corr'] = np.std([a['ph-caus-corr'] for a in last], ddof = 1)
    temp['ph-caus-ARI'] = np.std([a['ph-caus-ARI'] for a in last], ddof = 1)
    temp['ph-amp-caus-corr'] = np.std([a['ph-amp-caus-corr'] for a in last], ddof = 1)
    temp['ph-amp-caus-ARI'] = np.std([a['ph-amp-caus-ARI'] for a in last], ddof = 1)
    results.append(temp)

print total, modelscount

with open('CMIP5-interactions.csv', 'w') as f:
    writer = csv.writer(f, delimiter = ',')
    writer.writerow(['model', 'L1', 'L2', 'corr', 'synch corr', 'synch ARI', 'ph caus corr', 'ph caus ARI', 'ph amp corr', 'ph amp ARI'])
    for model, res in zip(models, results):
        writer.writerow([model, res['wvltL1'], res['wvltL2'], res['corr'], res['synch-corr'], res['synch-ARI'], res['ph-caus-corr'], res['ph-caus-ARI'], res['ph-amp-caus-corr'], res['ph-amp-caus-ARI']])