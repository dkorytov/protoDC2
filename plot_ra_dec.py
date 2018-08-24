#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr 
import h5py 
import dtk
import sys
from numpy.random import choice

param = dtk.Param(sys.argv[1])
lc_loc = param.get_string('lc_loc')
gltcs_steps = param.get_int_list('gltcs_steps')


color = 'b'

# plt.figure()
# for step in gltcs_steps:
#     print step
#     ra = np.fromfile(lc_loc.replace('${step}',str(step)).replace('${var_name}','theta'),dtype='f4')/3600
#     dec = np.fromfile(lc_loc.replace('${step}',str(step)).replace('${var_name}','phi'),dtype='f4')/3600
#     plt.plot(ra,dec,',b',alpha=0.005,)
# plt.ylabel('RA [deg]')
# plt.xlabel('Dec [deg]')
# plt.tight_layout()




plt.figure(figsize=(9,3))
y_bins = np.linspace(0,230,200)
x_bins = np.linspace(0,2600,2000)
H,_,_ = np.histogram2d([],[],bins=(x_bins,y_bins))
for step in gltcs_steps:
    print step

    x = np.fromfile(lc_loc.replace('${step}',str(step)).replace('${var_name}','x'),dtype='f4')
    y = np.fromfile(lc_loc.replace('${step}',str(step)).replace('${var_name}','y'),dtype='f4')
    # H2,_,_ = np.histogram2d(x,y,bins=(x_bins,y_bins))
    # H = H+H2
    slct = np.arange(x.size,dtype=int)
    slct = choice(slct,x.size/4,replace=False)
    plt.plot(x[slct],y[slct],'b,',alpha=0.002) 
plt.ylabel('y [Comoving Mpc/h]')
plt.xlabel('x [Comoving Mpc/h]')
plt.grid()
plt.xlim([0,2600])
plt.tight_layout()
    



# plt.figure()
# plt.pcolor(x_bins,y_bins,H.T,cmap='PuBu',norm=clr.LogNorm())
# plt.ylabel('y [Comoving Mpc/h]')
# plt.xlabel('x [Comoving Mpc/h]')
# plt.tight_layout()

dtk.save_figs('figs/'+__file__+'/'+param.file+'/')
plt.show()
