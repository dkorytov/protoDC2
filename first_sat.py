#!/usr/bin/env python2.7


import numpy as np
import matplotlib as plt
import dtk
import os.path

file_loc = '/home/dkorytov/data/AlphaQ/core_catalog3/01_12_17.AlphaQ.${step}.coreproperties'
steps = np.arange(0,500,dtype = int)

for step in steps:
    file_loc2 = file_loc.replace('${step}',str(step))
    if(os.path.isfile(file_loc2)):
        infall_step = dtk.gio_read(file_loc2,'infall_step')
        slct = infall_step != step
        print "\t",step
        print "\t\t",np.sum(slct)

        
        
