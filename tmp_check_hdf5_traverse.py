#!/usr/bin/env python2.7

import numpy as np
import h5py

gltc_loc = "output/mock_full_nocut_dust_elg_shear_mod.hdf5"


with h5py.File(gltc_loc, 'r') as fh:
    hgroup = fh['galaxyProperties']
    hobjects = []
    hgroup.visit(hobjects.append)
    hdatasets = [hobject for hobject in hobjects if type(hgroup[hobject])==h5py.Dataset]
    for data in hdatasets:
        print data
    native_quantities = set(hdatasets)

print " ===="
for nat in native_quantities:
    print type(nat),nat
    if(nat =='SEDs'):
        print '\t\tproblem'
    if(nat =='LSST_filters'):
        raise
