#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import dtk

param = dtk.Param(sys.argv[1])

output=param.get_string("output")
output_mod = output.replace(".hdf5","_mod.hdf5")
output_lsst = output.replace(".hdf5","_lsst.hdf5")

hfile = h5py.File(output_mod,'r')
hfile_lsst = h5py.File(output_lsst,'w')
keys = [
    'totalLuminositiesStellar:LSST_g:observed:dustAtlas',
    'totalLuminositiesStellar:LSST_r:observed:dustAtlas',
    'totalLuminositiesStellar:LSST_i:observed:dustAtlas',
    'totalLuminositiesStellar:LSST_z:observed:dustAtlas',
    'totalLuminositiesStellar:LSST_y:observed:dustAtlas',]

names = ['/lsst_g','/lsst_r','/lsst_i','/lsst_z','/lsst_y']

for i in range(0,len(keys)):
    a =keys[i]
    print hfile['/galaxyProperties/'].keys()
    b=hfile['/galaxyProperties/LSST_filters/'][a].value
    hfile_lsst[names[i]]=b

keys = ['redshift','redshiftHubble']
names = ['redshift_observed','redshift_cosmological']
for key,names in zip(keys,names):
    hfile_lsst[names]=hfile['/galaxyProperties/'+key].value
