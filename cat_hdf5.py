#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from astropy.cosmology import WMAP7 as cosmo
import sys
import time
from scipy.interpolate import interp1d 
import dtk

pd.set_option('display.expand_frame_repr', False)
h=0.702
total_start = time.time()
param = dtk.Param(sys.argv[1])

lc_loc = param.get_string("lc_loc")
gltcs_loc = param.get_string("gltcs_loc")
shear_loc = param.get_string("shear_loc")
gio_loc   = param.get_string("gio_loc")
use_shear = param.get_bool("use_shear")
gltcs_file_num = param.get_int("gltcs_file_num")
gltcs_file_step = param.get_int("gltcs_file_step")
gltcs_steps     = param.get_int_list("gltcs_steps")
gltcs_internal  = param.get_int_list("gltcs_internal")
gltcs_str_z     = param.get_string_list("gltcs_str_z")
gltcs_file_list = param.get_string_list("gltcs_file_list")
nan_str_z       = param.get_string_list("nan_str_z")
use_mr_cut      = False;
mr_cut          = param.get_float("mr_cut")
tmp_to_disk     = param.get_bool("tmp_to_disk")
tmp_file_loc    = param.get_string("tmp_file_loc")
output          = param.get_string("output")
stepz = dtk.StepZ(200,0,500)
zs = np.linspace(0,1.5,1000)
z_to_dl = interp1d(zs,cosmo.luminosity_distance(zs))
step2nan_str_z = {}

def cat_hdf5(file_loc,nums,output):
    keys = h5py.File(file_loc.replace("${num}","0"),'r').keys()
    hfile= h5py.File(output,'w')
    i_key_max = len(keys)
    for i_key, key in enumerate(keys):
        data = []
        print "%.2f  %d/%d"%(float(i_key)/float(i_key_max),i_key,i_key_max)
        for i in range(0,nums):
            htmp = h5py.File(file_loc.replace("${num}",str(i)),'r')
            data.append(htmp[key][:])
        data_tot = np.concatenate(data)
        hfile[key]=data_tot
        hfile.flush()

tmp_file_num = len(gltcs_file_list)
print "cat from", tmp_file_loc," num files: ",tmp_file_num," to ",output
cat_hdf5(tmp_file_loc,tmp_file_num,output)
