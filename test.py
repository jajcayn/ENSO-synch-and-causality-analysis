import numpy as np
import matplotlib.pyplot as plt
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
from matplotlib.ticker import MultipleLocator, FuncFormatter
import cPickle
import sys
import matplotlib.gridspec as gridspec
sys.path.append('/home/nikola/Work/phd/multi-scale/src')
sys.path.append('/home/nikola/Work/phd/multi-scale/surrogates')
sys.path.append('/home/nikola/Work/phd/mutual_information')

from data_class import DataField
from surrogates import SurrogateField
from multiprocessing import Process, Queue
from time import sleep

CMIP5model = None

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

a = list(enso.get_seasonality(DETREND = False))
enso_sg = SurrogateField()

_, _, idx = enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = False)
enso_sg.copy_field(enso)
enso_sg.data = enso_sg.data[idx[0]:idx[1]]

enso.return_seasonality(a[0], a[1], None)

a[0] = a[0][idx[0]:idx[1]]
a[1] = a[1][idx[0]:idx[1]]
enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = True)

enso.center_data()
print 'data', enso.data[:5]

def _surrs(sg, a, jobq, resq):
    mean, var, _ = a
    while jobq.get() is not None:
        # np.random.seed()
        sg.construct_fourier_surrogates_spatial()
        n = np.random.uniform(0, 2*np.pi)
        ts = sg.get_surr()
        # sg.add_seasonality(mean, var, None)

        # sg.center_surr()

        sleep(1)
        
        resq.put((ts, n))

numsurr = 100
WRKRS = 20
comp = 0
jobq = Queue()
resq = Queue()
for i in range(numsurr):
    jobq.put(1)
for i in range(WRKRS):
    jobq.put(None)
wrkrs = [Process(target=_surrs, args = (enso_sg, a, jobq, resq)) for i in range(WRKRS)]
for w in wrkrs:
    w.start()
    
while comp < numsurr:
    ts, n = resq.get()
    print ts[:5]
    # print n
    comp += 1
    
for w in wrkrs:
    w.join()