#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas
import sys
import dtk


def add_attrs(hfile,attrs):
    for attr in attrs.keys():
        hfile.attrs[attr] = attrs[attr]


stepz = dtk.StepZ(200,0,500)
output = "output/mock_shear_full.hdf5"
print "loading data"
cat_pd = dtk.pandas_from_hdf5(output)
print "done loading"
steps = np.unique(cat_pd['step'])
attrs = {}
hfile = h5py.File(output,'r')
attrs_list = hfile.attrs.keys()
print attrs_list
for key in attrs_list:
    attrs[key] = hfile.attrs[key]
hfile = h5py.File("output/bystep/step_mock.hdf5",'w')
add_attrs(hfile,attrs)
for step in steps:
    print step
    slct = cat_pd['step']==step
    output_step = "output/bystep/"+str(step)+"_mock.hdf5"
    step_cat_pd = cat_pd[slct]
    dtk.pandas_to_hdf5(step_cat_pd,output_step)
    hfile_tmp = h5py.File(output_step,'a')
    add_attrs(hfile_tmp,attrs)
    for key in step_cat_pd.keys():
        hfile[""+str(step)+"/"+key] = step_cat_pd[key]
    hfile[""+str(step)+"/"].attrs['redshift'] = stepz.get_z(step)
    add_attrs(hfile,attrs)
hfile.close()
