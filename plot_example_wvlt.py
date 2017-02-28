import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from datetime import date, datetime
# from matplotlib.ticker import MultipleLocator, FuncFormatter
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
sys.path.append('/home/nikola/Work/multi-scale')
from src.data_class import DataField, load_enso_index, load_station_data


FROM = 1970
TO = 2000

enso = load_enso_index("nino34raw.txt", '3.4', date(1900, 1, 1), date(2014, 1, 1), anom = True)
ndx = enso.select_date(date(FROM-1, 1, 1), date(TO+1, 1, 1), apply_to_data = True)


plt.figure(figsize=(15,10))
# plt.figure(figsize=(15,7))
plt.subplot(2, 1, 1)
enso.wavelet(2, 'y', cut = 1)
plt.xlim([0, enso.phase.shape[0]])
plt.ylim([-4,4])
plt.plot(enso.phase, linewidth = 2)
plt.plot(enso.amplitude, linewidth = 2)
plt.plot(enso.data[12:-12], color = '#242632', linewidth = 2.5)
x = np.arange(enso.find_date_ndx(date(1982, 9, 1)), enso.find_date_ndx(date(1983, 5, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1997, 7, 1)), enso.find_date_ndx(date(1998, 4, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1972, 8, 1)), enso.find_date_ndx(date(1973, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1973, 9, 1)), enso.find_date_ndx(date(1974, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1975, 9, 1)), enso.find_date_ndx(date(1976, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1988, 8, 1)), enso.find_date_ndx(date(1989, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
plt.ylabel("$A_{2}$, $\phi_{2}$", size = 28)
plt.xticks(np.arange(0, enso.data.shape[0]+1, ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10), rotation = 30)
plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(19) 
plt.gca().spines['top'].set_visible(False)
# plt.grid()
plt.gca().spines['bottom'].set_visible(False)    
plt.gca().spines['right'].set_visible(False) 

plt.subplot(2, 1, 2)
enso.wavelet(2, 'y', cut = 1)
plt.plot(enso.amplitude * np.cos(enso.phase), linewidth = 2)
enso.wavelet(1, 'y', cut = 1)
plt.plot(enso.amplitude * np.cos(enso.phase), linewidth = 2)
plt.plot(enso.data[12:-12], color = '#242632', linewidth = 2.5)
plt.xlim([0, enso.phase.shape[0]])
plt.ylim([-4, 4])
x = np.arange(enso.find_date_ndx(date(1982, 9, 1)), enso.find_date_ndx(date(1983, 5, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1997, 7, 1)), enso.find_date_ndx(date(1998, 4, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1972, 8, 1)), enso.find_date_ndx(date(1973, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#F83F67", edgecolor = "#F83F67", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1973, 9, 1)), enso.find_date_ndx(date(1974, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1975, 9, 1)), enso.find_date_ndx(date(1976, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
x = np.arange(enso.find_date_ndx(date(1988, 8, 1)), enso.find_date_ndx(date(1989, 3, 1)), 1)
plt.fill_between(x, -4, 4, facecolor = "#30AEDF", edgecolor = "#30AEDF", alpha = 0.5)
plt.ylabel("$A_{2}\cos{\phi_{2}}$, $A_{1}\cos{\phi_{1}}$", size = 28)
plt.xticks(np.arange(0, enso.data.shape[0]+1, ((TO-FROM)/10)*12), np.arange(FROM, TO+1, (TO-FROM)/10), rotation = 30)
plt.gca().xaxis.set_minor_locator(MultipleLocator(12))
for tick in plt.gca().xaxis.get_major_ticks():
    tick.label.set_fontsize(19)
for tick in plt.gca().yaxis.get_major_ticks():
    tick.label.set_fontsize(19) 
plt.gca().spines['top'].set_visible(False)
# plt.grid()
plt.gca().spines['bottom'].set_visible(False)    
plt.gca().spines['right'].set_visible(False) 


# plt.show()
plt.savefig("wvlt_example_double.eps", bbox_inches = "tight")

