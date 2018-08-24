#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import h5py 
import dtk

hfile = h5py.File('output/mock_full_nocut_dust_elg_shear3_mod.hdf5','r')

conv = hfile['galaxyProperties/convergence'].value
mag = hfile['galaxyProperties/magnification'].value
shear1 = hfile['galaxyProperties/shear1'].value
shear2 = hfile['galaxyProperties/shear2'].value
ra  = hfile['galaxyProperties/ra'].value
dec  = hfile['galaxyProperties/dec'].value
conv_new = 1.0-np.sqrt(1.0/mag + shear1**2 + shear2**2)
print conv_new
xbins = np.linspace(-.05,.1,100)
plt.figure()
h,xbins = np.histogram(conv, bins = xbins)
plt.plot(dtk.bins_avg(xbins),h,'-',label='old')
h,xbins = np.histogram(conv_new, bins =xbins)
plt.plot(dtk.bins_avg(xbins),h,'-',label='new')
plt.legend(loc='best')

print np.min(mag), np.max(mag)


plt.figure()
slct  = mag < 0.5
xbins = np.linspace(-.1,0.5,100)
h,xbins = np.histogram(mag[slct], bins = xbins)
plt.plot(dtk.bins_avg(xbins),h)

plt.figure()
plt.plot(ra[slct],dec[slct],'.',alpha=0.3)
plt.show()

