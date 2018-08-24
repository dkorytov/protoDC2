#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import dtk
import sys
import myutil
import h5py

param1 = dtk.Param(sys.argv[1])
param2 = dtk.Param(sys.argv[2])


output1 = param1.get_string('output')
output_mod1 = output1.replace(".hdf5","_mod.hdf5")
hfile1 =  h5py.File(output_mod1,'r')

output2 = myutil.get_output_loc(param2)
output_mod2 = output2.replace(".hdf5","_mod.hdf5")
hfile2 =  h5py.File(output_mod2,'r')

match_list = []
mismatch_list = []
mismatch_val  = []
only_1_list = []
only_2_list = []

def check_group(hgroup1,hgroup2):
    hobjects1 = []
    hobjects2 = []
    hgroup1.visit(hobjects1.append)
    hgroup2.visit(hobjects2.append)
    keys1 = [hobject for hobject in hobjects1 if type(hgroup1[hobject]) == h5py.Dataset]
    keys2 = [hobject for hobject in hobjects2 if type(hgroup2[hobject]) == h5py.Dataset]
    for key in keys1:
        if key not in keys2:
            print "not in 2: ",key
            only_1_list.append(key)
    for key in keys2:
        if key not in keys1:
            print "not in 1: ",key
            only_2_list.append(key)

    for key in keys1:
        print key
        if key in keys2:
            data1 = hgroup1[key].value
            data2 = hgroup2[key].value
            diff = data1-data2
            if(np.sum(diff!=0)>0):
                mismatch_list.append(key)
                mismatch_val.append(diff)
            else:
                match_list.append(key)

check_group(hfile1['galaxyProperties'],hfile2['galaxyProperties'])

print "only1: ",len(only_1_list)
print "only2: ",len(only_2_list)
print "match_list size: ",len(match_list)
print "mismatch_list size: ",len(mismatch_list)

            
