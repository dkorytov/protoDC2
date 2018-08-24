#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import dtk
import h5py
import myutil

param = dtk.Param(sys.argv[1])
#output=param.get_string('output')
sim_step = param.get_int_list('gltcs_steps')
glctcs_step = param.get_int_list('gltcs_internal')
output = myutil.get_output_loc(param)
output_mod = output.replace('.hdf5','_mod.hdf5')
glctcs_files = param.get_string_list('gltcs_file_list')
glctcs_count = np.zeros(len(glctcs_files))
# bad_gal_loc = sys.argv[2]

# bad_gal = np.loadtxt(bad_gal_loc,dtype={'names':('id','av','bv','mag_g'),
#                                         'formats':('i4','f8','f8','f8')})
# bad_gal_id = np.concatenate((bad_gal['id'], [21648285,17869188 ]))
bad_gal_id = [14102637, 1224979446537128410, 1418634088780137682, 873698490918637197]
s2g = {}
for ss,gs in zip(sim_step,glctcs_step):
    s2g[ss]=gs

hfile = h5py.File(output_mod,'r')
galID = hfile['galaxyProperties/galaxyID'].value
nodeIndex = hfile['galaxyProperties/infallIndex'].value
step = hfile['galaxyProperties/step'].value
#id_list = [21648285,17869188 ]
#id_list = [17869188 ]

indx = dtk.search_sorted(galID,bad_gal_id)
#indx = dtk.search_sorted(galID,id_list)
print step[indx][:10]
print galID[indx][:10]
# print bad_gal['id'][:10]

print step[indx]
bad_nodeIndex = nodeIndex[indx]
bad_step = step[indx]
nodeIndexs_to_find = nodeIndex[indx]
max_bad = 0
for i,glctcs_file in enumerate(glctcs_files):
    print i,'/',len(glctcs_files), glctcs_file
    glctcs_nodeIndex = h5py.File(glctcs_file,'r')['Outputs/Output6/nodeData/nodeIndex'].value
    srt = np.argsort(glctcs_nodeIndex)
    indexs = dtk.search_sorted(glctcs_nodeIndex,bad_nodeIndex,sorter=srt)
    found_num = np.sum(indexs != -1)
    if(found_num > max_bad):
        max_bad = found_num
    print '\t found bad: ', found_num
    # if(21648285 in bad_gal_id[indexs!=-1]):
    #     print "this file: "
    #     print glctcs_file
    #     slct = indexs != -1
    #     for bdgalid, bindx,bdstep in zip(bad_gal_id[slct],bad_nodeIndex[slct],bad_step[slct]):
    #         print bdgalid,bindx,bdstep,s2g[bdstep]
    #     break

print max_bad

