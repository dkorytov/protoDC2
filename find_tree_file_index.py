#!/usr/bin/env python2.7

import numpy as np
import h5py 
import dtk
import sys

param = dtk.Param(sys.argv[1])

merger_tree_loc = param.get_string("merger_tree_loc")
merger_tree_num = param.get_int("merger_tree_num")
search_nodeIndexs = param.get_int64_list("nodeIndexs")


for file_num in range(64,merger_tree_num):
    file_name = merger_tree_loc.replace("${subfile}",str(file_num))
    hfile = h5py.File(file_name);
    firstNodes = hfile['/forestIndex/firstNode'].value
    numNodes   = hfile['/forestIndex/numberOfNodes'].value
    #forestIndexs = hfile['/forestIndex/forestIndex'].value
    nodeIndexs = hfile['/forestHalos/nodeIndex'].value
    
    indx = np.in1d(nodeIndexs,search_nodeIndexs).nonzero()[0]
    if(indx.size > 0):
        firstNode_num = np.searchsorted(firstNodes,indx)-1
        for j in range(0,indx.size):
            for i in range(0,search_nodeIndexs.size):
                if(search_nodeIndexs[i]==nodeIndexs[indx[j]]):
                    print "nodeIndex: ", search_nodeIndexs[i]
                    print '\tfile num: ',file_num
                    print '\ttree num: ',firstNode_num[j]
                    #testing
                    found = False
                    for k in range(0,numNodes[firstNode_num[j]]):
                        if(nodeIndexs[firstNodes[firstNode_num[j]]+k]==search_nodeIndexs[i]):
                            found = True
                    if(not found):
                        print "Not found?! wtf?!"
                        
    # for i in range(0,firstNodes.size):
    #     for j in range(0,numNodes[i]):
    #         nindx =nodeIndexs[firstNodes[i]+j] 
    #         for k in range(0,search_nodeIndexs.size):
    #             if( nindx == search_nodeIndexs[k]):
    #                 print "node index: ", search_nodeIndexs[k]
    #                 print "\n\ttree num:", i
                                   
