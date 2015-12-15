import cPickle
import numpy as np
import matploltlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt




fname = "zCMI_PROdampedSST_1000FTsurrs_16384_ENSOperiod=4.75_lambda0=2.bin"
NUM_SURR = 1000

with open(fname, "rb") as f:
    data = cPickle.load(f)


print data.keys()

MI = data['zMI']
CMI1 = data['zCMI1']
CMI2 = data['zCMI2']



SPAN = [0.5, 7.5]
y = 12
scales = np.arange(SPAN[0]*y, SPAN[1]*y+1, 1)


plt.figure(figsize = (15,10))
p = plt.subplot(311)
p.tick_params(axis='both', which='major', labelsize = 17)
p1, = plt.plot(scales, MI, color = "#0059C7", linewidth = 2)
plt.xlim([scales[0], scales[-1]])
plt.ylim([0, 10])
if NUM_SURR > 0:
    plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
for s in scales[6::12]:
    plt.axvline(x = s, color = "#cccccc", linewidth = 0.6)
plt.legend([p1], ['MI // 4bins -- phaseAnn vs. phaseOthers'])
plt.xticks(scales[6::12], [int(s) for s in scales[6::12] / 12])

p = plt.subplot(312)
p.tick_params(axis='both', which='major', labelsize = 17)
p1, = plt.plot(scales, CMI2, color = "#0059C7", linewidth = 2)
plt.xlim([scales[0], scales[-1]])
plt.ylabel("information [bits]", size = 25)
plt.ylim([-2, 14])
plt.legend([p1], ['CMI // 4 bins -- phaseAnn -> phaseOthers'])
plt.xticks(scales[6::12], [int(s) for s in scales[6::12] / 12])
for s in scales[6::12]:
    plt.axvline(x = s, color = "#cccccc", linewidth = 0.6)
if NUM_SURR > 0:
    plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
    plt.ylabel("z-score", size = 25)

p = plt.subplot(313)
p.tick_params(axis='both', which='major', labelsize = 17)
p1, = plt.plot(scales, CMI1, color = "#0059C7", linewidth = 2)
plt.xlim([scales[0], scales[-1]])
plt.ylim([-2, 20])
plt.legend([p1], ['CMI // 4bins -- phaseOthers -> phaseAnn'])
plt.xticks(scales[6::12], [int(s) for s in scales[6::12] / 12])
for s in scales[6::12]:
    plt.axvline(x = s, color = "#cccccc", linewidth = 0.6)
if NUM_SURR > 0:
    plt.axhline(y = 2., color = "#910D3E", linewidth = 0.75)
plt.xlabel("period [years]", size = 20)

# plt.suptitle("NINO3.4 SST // z-score against %d FT surrogates" % (NUM_SURR), size = 25)
plt.suptitle("PRO damped model // length = 16384 // ENSO period = 4.75 //\n $\lambda_{0}$ = 2yr$^{-1}$ // z-score against %d FT surrogates" % (NUM_SURR), size = 30)
# plt.suptitle("PRO model -- neutral SST // z-score against %d FT surrogates" % (NUM_SURR), size = 25)
# plt.show()
plt.savefig('test.png')