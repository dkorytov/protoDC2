#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import dtk


param = dtk.Param(sys.argv[1])

output = param.get_string("output")
output_mod = output.replace('.hdf5','_mod.hdf5')

print output
hfile = h5py.File(output_mod,'r')
dec =  hfile['galaxyProperties/ra'][:]
ra = hfile['galaxyProperties/dec'][:]
m0 =  hfile['galaxyProperties/magnification'][:]
k0 =  hfile['galaxyProperties/convergence'][:]
s1 =  hfile['galaxyProperties/shear1'][:]
s2 =  hfile['galaxyProperties/shear2'][:]
ss = np.sqrt(s1**2 + s2**2)
mr = hfile['galaxyProperties/SDSS_filters/magnitude:SDSS_r:rest'][:]
redshift = hfile['galaxyProperties/redshift'][:]
order=np.argsort(ss)
order2 = np.argsort(m0)
slct = mr<-21.5
print "mr cut: ", float(slct.sum())/float(slct.size)
plt.figure()
plt.plot(ra[slct]/60/60,dec[slct]/60/60,',',alpha=0.1)
plt.ylabel('dec')
plt.xlabel('ra')
plt.figure()
ax = plt.subplot(111,projection='polar')
print np.min(redshift)
plt.plot(dec[slct]/60/60/180*np.pi,redshift[slct],',',alpha=0.005)
rs = np.linspace(0.2,1.2,7)
rs_label = ["z=%.1f"%(r) for r in rs]
plt.plot(np.pi,1.4,'x')
ax.set_rgrids(rs,labels=rs_label,angle=7.5)
plt.xlabel('angle')

plt.figure()
ax = plt.subplot(111,projection='polar')
ax.scatter(dec[order2][slct]/60/60/180*np.pi,redshift[order2][slct],c=np.log(m0[order2][slct]-1),s=5,edgecolors='none',alpha=0.005)
ax.set_xlabel('angle')
plt.plot(np.pi,1.3,'x')
rs = np.linspace(0.2,1.2,7)
rs_label = ["z=%.1f"%(r) for r in rs]
ax.set_rgrids(rs,labels=rs_label,angle=7.5)



plt.figure()
plt.scatter(ra[order][slct]/60/60,dec[order][slct]/60/60,c=ss[order][slct],s=5,edgecolors='none',alpha=0.5)
plt.ylabel('Dec')
plt.xlabel('RA')
cb = plt.colorbar()
cb.set_label('shear magnitude')
dtk.save_figs('figs/'+param.file+'/'+__file__+'/')
plt.show()
