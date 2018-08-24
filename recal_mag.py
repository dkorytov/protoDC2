#!/usr/bin/env python2.7



#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import h5py 
import dtk
import sys

from astropy.cosmology import WMAP7 as cosmo
import time
from scipy.interpolate import interp1d 


h=0.702
stepz = dtk.StepZ(200,0,500)
zs = np.linspace(0,1.5,1000)
z_to_dl = interp1d(zs,cosmo.luminosity_distance(zs))
step2nan_str_z = {}

param = dtk.Param(sys.argv[1])
output = param.get_string('output')
output_mod = output.replace('.hdf5','_mod.hdf5')

hfile = h5py.File(output_mod,'a')

redshift = hfile['galaxyProperties/redshiftHubble'].value
dl = z_to_dl(redshift)
print "done getting the luminosity distance..."
adjust = -2.5*np.log10(1+redshift)+5*np.log10(dl)+25.0



def recal_obs(hgroup):
    keys = hgroup.keys()
    new_val = []
    old_val = []
    for key in keys:
        if ('diskLuminositiesStellar' in key and ('observed' in key)):
            key_base = key.replace('diskLuminositiesStellar','')
            key_mag = 'magnitude'+key_base
            lum_d = hgroup['diskLuminositiesStellar'+key_base].value
            lum_s = hgroup['spheroidLuminositiesStellar'+key_base].value
            new_val = adjust - 2.5*np.log10(lum_d+lum_s)
            old_val = hgroup[key_mag].value
            hgroup[key_mag][...]=new_val
            print key_base
            for i in range(0,100):
                print '%.2f -> %.2f' %(old_val[i],new_val[i])

recal_obs(hfile['galaxyProperties/SDSS_filters'])
recal_obs(hfile['galaxyProperties/LSST_filters'])
