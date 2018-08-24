#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import dtk
import h5py
import myutil

param = dtk.Param(sys.argv[1])
glctcs_files = param.get_string_list('gltcs_file_list')
glctcs_count = np.zeros(len(glctcs_files))

node_list = []
for i in range(2,len(sys.argv)):
    node_list.append(int(sys.argv[i]))




nodeIndexs_to_find = np.array(node_list)
print nodeIndexs_to_find
max_bad = 0
for i,glctcs_file in enumerate(glctcs_files):
    print '\n', i,'/',len(glctcs_files), glctcs_file
    for j in range(1,33):
        print "\t",j,
        glctcs_nodeIndex = h5py.File(glctcs_file,'r')['Outputs/Output%d/nodeData/nodeIndex'%(j)].value

        srt = np.argsort(glctcs_nodeIndex)
        indexs = dtk.search_sorted(glctcs_nodeIndex,nodeIndexs_to_find,sorter=srt)
        found_num = np.sum(indexs != -1)
        if(found_num>0):
            print found_num
            raise

print max_bad

