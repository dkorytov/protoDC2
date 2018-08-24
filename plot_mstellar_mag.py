#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import h5py
import dtk
import sys

param = dtk.Param(sys.argv[1])
output= param.get_string("output")
output_mod = output.replace('.hdf5','_mod.hdf5')
hfile = h5py.File(output_mod,'r')
print("loading stellar_mass")
stellar_mass = np.log10(hfile['galaxyProperties/totalMassStellar'].value)
print("loading mag_r")
mag_r = hfile['galaxyProperties/SDSS_filters/magnitude:SDSS_r:rest'].value

H,x_bins,y_bins = np.histogram2d(stellar_mass,mag_r,bins=(100,100))

plt.figure()
plt.pcolor(x_bins,y_bins,H.T,cmap='PuBu',norm=clr.LogNorm())
cb = plt.colorbar()
plt.ylabel('SDSS r-band abs. Mag.')
plt.xlabel('Log10(Stellar Mass [Msun/h])')
cb.set_label('Freq')
plt.tight_layout()
plt.grid()

plt.figure()
plt.pcolor(x_bins,y_bins,H.T,cmap='PuBu')
cb = plt.colorbar()
plt.ylabel('SDSS r-band abs. Mag.')
plt.xlabel('Log10(Stellar Mass [Msun/h])')
cb.set_label('Freq')
plt.tight_layout()
plt.grid()

H,x_bins = np.histogram(stellar_mass,bins=100)
x_bins_avg = dtk.bins_avg(x_bins)
plt.figure()
plt.plot(x_bins_avg,H)
plt.ylabel('Freq')
plt.xlabel('Log10(Stellar Mass [Msun/h])')
plt.grid()

H,x_bins = np.histogram(mag_r,bins=100)
x_bins_avg = dtk.bins_avg(x_bins)
plt.figure()
plt.plot(x_bins_avg,H)
plt.ylabel('Freq')
plt.xlabel('SDSS r-band abs. Mag.')

plt.grid()
# plt.figure()
# plt.pcolor(x_bins,y_bins,H.T,cmap='viridis',norm=clr.LogNorm())
# plt.colorbar()
dtk.save_figs('figs/'+param.file+'/'+__file__+'/')

plt.show()



