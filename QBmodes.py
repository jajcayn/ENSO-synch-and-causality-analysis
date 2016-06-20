import numpy as np
import sys
sys.path.append('/Users/nikola/work-ui/multi-scale')
from src.data_class import DataField
import matplotlib.pyplot as plt
import scipy.signal as ss
from mpl_toolkits.basemap import Basemap
from mpl_toolkits import basemap
import matplotlib.pyplot as plt

def plot(field, lats, lons, u = None, v = None, tit = None, fname = None):
    plt.figure(figsize=(20,15))
    lat_ndx = np.argsort(lats)
    lats = lats[lat_ndx]
    m = Basemap(projection = 'merc',
                    llcrnrlat = lats[0], urcrnrlat = lats[-1],
                    llcrnrlon = lons[0], urcrnrlon = lons[-1],
                    resolution = 'c')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,360.,60.))
    x, y = m(*np.meshgrid(lons, lats))
    m.contourf(x, y, field[lat_ndx, :])
    if tit is not None:
        plt.title(tit, size = 25)
    if u is not None and v is not None:
        m.quiver(x, y, u, v)
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = 'tight')
    plt.close()


def pcs_spectra(g, num = 10):
    eofs, pcs, var = g.pca_components(n_comps = num)
    var *= 100
    for i in range(num):
        f, pxx = ss.welch(pcs[i, :], fs = 1./(3*2.628e+6), scaling = 'spectrum')
        f *= 3.154e+7 # in 1/years
        plt.figure()
        plt.semilogx(f, np.sqrt(pxx), color = 'k', linewidth = 2)
        plt.xlim([0, 6])
        plt.xlabel("1/years", size = 20)
        plt.ylabel("power", size = 20)
        plt.fill_between(np.arange(1/2.5, 1/1.5, 0.01), 0, plt.yticks()[0][-1], facecolor = "grey", alpha = 0.6)
        plt.title("%s - PC %d spectra (%.3f %%)" % (g.var_name, i+1, var[i]))
        plt.show()
        # plt.savefig("plots/QB/%s-PCspectra%d.png" % (g.var_name, i+1))
        plt.close()
        # plot(eofs[i, ...], g.lats, g.lons, tit = "%s - EOF%d (%.3f%%)" % (g.var_name, i+1, var[i]), fname = "plots/QB/%s-EOF%d.png" % (g.var_name, i+1))


folder = '/Users/nikola/work-ui/data/'

# sst = DataField(data_folder = folder)
# sst.load(filename = '20CR.sst.mon.nc', variable_name = 'sst', dataset = 'ERA')
# sst.select_lat_lon([-30, 30], [100, 300])
# sst.anomalise()

# slp = DataField(data_folder = folder)
# slp.load(filename = '20CR.slp.mon.nc', variable_name = 'msl', dataset = 'ERA')
# slp.select_lat_lon([-30, 30], [100, 300])
# slp.anomalise()

# u10 = DataField(data_folder = folder)
# u10.load(filename = '20CR.wind10m.mon.nc', variable_name = 'u10', dataset = 'ERA')
# v10 = DataField(data_folder = folder)
# v10.load(filename = '20CR.wind10m.mon.nc', variable_name = 'v10', dataset = 'ERA')
# u10.select_lat_lon([-30, 30], [100, 300])
# v10.select_lat_lon([-30, 30], [100, 300])
# u10.anomalise()
# v10.anomalise()

# pcs_spectra(v10)

td = DataField(data_folder = folder)
td.load(filename = 'heat_content_anomaly_0-700_seasonal.nc', variable_name = 'h18_hc', dataset = 'NCEP')
td.select_level(0)
td.select_lat_lon([-30, 30], [150, 250])
# td.anomalise()
pcs_spectra(td)
