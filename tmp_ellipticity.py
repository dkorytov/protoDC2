#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import dtk
import sys
import h5py
from numpy.random import rand
from numpy.random import seed
from scipy.interpolate import interp1d

param = dtk.Param(sys.argv[1])
output = param.get_string("output")
hfile = h5py.File(output,'r')


def gal_zoo_dist(x):
    val = np.zeros_like(x)
    slct = x<0.2
    val[slct] = 0

    slct = (0.2<=x) & (x<0.4)
    val[slct] = (x[slct]-0.2)*5

    slct = (0.4<=x) & (x<0.8)
    val[slct] = 1.0

    slct = (0.8<=x) & (x<1.0)
    val[slct] = (1-x[slct])*5.0

    slct = 1.0<=x
    val[slct] = 0
    return val

def gal_zoo_dist(x):
    val = np.zeros_like(x)
    a = 2
    slct = x<0.2
    val[slct] = 0

    slct = (0.1<=x) & (x<0.7)
    val[slct] = np.tanh((x[slct]-.3)*np.pi*a) - np.tanh((-0.2)*np.pi*a)

    slct = (0.7<=x) & (x<1.0)
    val[slct] = np.tanh(-(x[slct]-.95)*np.pi*6.) - np.tanh((-0.2)*np.pi*a) -(np.tanh(-(0.7-0.95)*np.pi*6)-np.tanh(0.4*np.pi*a))

    slct = 1.0<=x
    val[slct] = 0
    return val

plt.figure()
x=np.linspace(0,1,1000)
plt.plot(x,gal_zoo_dist(x))
#plt.plot(x,gal_zoo_dist2(x))

inclination = hfile['inclination'].value/180*np.pi
size = inclination.size
pos_angle = np.array(rand(size)*np.pi)
q = dtk.clipped_gaussian(0.8, 0.2, size, max_val = 1.0, min_val=0.0)

print np.min(q),np.max(q)

spheroid_ellip = (1.0-q)/(1.0+q)
spheroid_ellip1 = np.cos(2.0*pos_angle)*spheroid_ellip
spheroid_ellip2 = np.sin(2.0*pos_angle)*spheroid_ellip
plt.figure()
h,xbins = np.histogram(q,100)
plt.plot(dtk.bins_avg(xbins),h,label='axis ratio')
h,xbins=np.histogram(spheroid_ellip,100)
plt.plot(dtk.bins_avg(xbins),h,label = 'ellipticity')
plt.legend(loc='best')


plt.figure()
h,xbins = np.histogram(spheroid_ellip1,100)
plt.plot(dtk.bins_avg(xbins),h,label='ellip1')
h,xbins = np.histogram(spheroid_ellip2,100)
plt.plot(dtk.bins_avg(xbins),h,label='ellip2')
plt.legend(loc='best')

plt.figure()
H,xbins,ybins=np.histogram2d(spheroid_ellip1,spheroid_ellip2,bins=(100,100))
plt.pcolor(xbins,ybins,H,cmap='PuBu',norm=clr.LogNorm())




# h,xbins = np.histogram(inclination,100)
# xbins_avg = dtk.bins_avg(xbins)
# plt.figure()
# plt.plot(xbins_avg,h)

dist,lim = dtk.make_distribution(-inclination)




# plt.figure()
# plt.plot(xbins_avg,gal_zoo_dist(xbins_avg))



resamp = dtk.resample_distribution(dist,gal_zoo_dist,lim,[0.0,1.0])


q_disk = resamp(-inclination)
disk_ellip = (1.0-q_disk)/(1.0+q_disk)
disk_ellip1 = np.cos(2.0*pos_angle)*disk_ellip
disk_ellip2 = np.sin(2.0*pos_angle)*disk_ellip
plt.figure()
h,xbins = np.histogram(q_disk,100)
plt.plot(dtk.bins_avg(xbins),h,label='disk axis ratio')
h,xbins = np.histogram(disk_ellip,100)
plt.plot(dtk.bins_avg(xbins),h,label='disk ellip')
plt.legend(loc='best')

plt.figure()
h,xbins = np.histogram(disk_ellip1,100)
plt.plot(dtk.bins_avg(xbins),h,label='disk ellip1')
h,xbins = np.histogram(disk_ellip2,100)
plt.plot(dtk.bins_avg(xbins),h,label='disk ellip2')
plt.legend(loc='best')




plt.show()


