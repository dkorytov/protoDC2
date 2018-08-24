#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as clr
import dtk
import h5py
import sys

param = dtk.Param(sys.argv[1])
output = param.get_string("output")
output_mod = output.replace(".hdf5","_mod.hdf5")
hfile =h5py.File(output_mod,'r')
keys = ['hostIndex','hostHaloMass','galaxyID','isCentral','step','redshiftHubble','redshift','x','y','z','vx','vy','vz']
gal_dict= {}
for key in keys:
    print key
    gal_dict[key] = hfile['galaxyProperties/'+key].value

print 'srting...'
srt = np.argsort(gal_dict['hostIndex'])
hfile_out = h5py.File(output.replace(".hdf5","_lc.hdf5"),'w');print "a"
for key in keys:
    print key
    hfile_out['galaxiesProperties/'+key]=gal_dict[key][srt]

slct = gal_dict['isCentral']==1
for key in keys:
    print key
    if key != 'isCentral':
        hfile_out['halos/'+key]=gal_dict[key][slct]
