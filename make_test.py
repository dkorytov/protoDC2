#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import sys
import h5py
import dtk


param = dtk.Param(sys.argv[1])
output = param.get_string("output")
output_mod = output.replace(".hdf5","_mod.hdf5")
output_test = output.replace(".hdf5","_micro.hdf5")

hfile = h5py.File(output_mod,'r')
hfile_test = h5py.File(output_test,'w')
ra = hfile['galaxyProperties/ra'].value
dec = hfile['galaxyProperties/dec'].value
#ra_len = hfile['galaxyProperties/ra_lensed'].value
#dec_len = hfile['galaxyProperties/dec_lensed'].value


print np.min(ra),np.max(ra)
print np.min(dec),np.max(dec)

# print np.min(ra_len),np.max(ra_len)
# print np.min(dec_len),np.max(dec_len)

slct = (-.1<ra) & (ra<.1) & (-.1<dec) & (dec<.1)
print np.sum(slct)

hfile_test.create_group('/galaxyProperties')
hgroup = hfile['galaxyProperties']
hobjects = []
#get all the names of objects in this tree
hgroup.visit(hobjects.append)
#filter out the group objects and keep the dataste objects
hdatasets = [hobject for hobject in hobjects if type(hgroup[hobject])==h5py.Dataset]
nat_quants = set(hdatasets)

print nat_quants
for nat_quant in nat_quants:
    print 'working on ', nat_quant
    hdataset = hfile['galaxyProperties/'+nat_quant]
    data = hdataset.value[slct]
    hfile_test['galaxyProperties/'+nat_quant] = data
    for attr in hdataset.attrs:
        hfile_test['galaxyProperties/'+nat_quant].attrs[attr] = hdataset.attrs[attr]

#copy meta data
hfile_test.copy(hfile['metaData'],'/metaData')
print hfile_test['metaData'].keys()
hfile_test['metaData/skyArea'][...]=1.0
#testing
ra1 = hfile_test['galaxyProperties/ra'].value
dec1 = hfile_test['galaxyProperties/ra'].value
print np.min(ra1),np.max(ra1)
print np.min(dec1),np.max(dec1)
