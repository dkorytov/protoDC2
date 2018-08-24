#!/usr/bin/env python2.7

import numpy as np
from astropy.cosmology import WMAP7 as cosmo
from scipy.interpolate import interp1d 
from scipy.integrate import cumtrapz
from numpy.random import rand
import matplotlib.pyplot as plt
import dtk
import sys
import pandas

output = "output/mock_shear_full.hdf5"
print "loading data"
cat=dtk.pandas_from_hdf5(output)
for i in range(0,10):
    print(cat['ra'][i]/60/60,cat['dec'][i]/60/60)

print "done loading"
def bt_func(disk,sphere):
    return sphere/(disk+sphere)

def sersic(r,n):
    return np.exp(-r**1.0/n)

def sersic_r(r,n):
    return r*np.exp(-r**1.0/n)

def sersic_inter(r,n):
    result = cumtrapz(sersic_r(r,n),x=r,initial=0.0)
    result = result/result[-1]
    return result

def sersic_half_scale_ratio(n):
    #half light radius is radius that contains 0.5 of the total light
    #scale radius is the radias when the profile intensity drops by 1/e. In this function, r_scale =1.0
    r = np.linspace(0,10,1000) 
    tot_light = sersic_inter(r,n)
    r_half = r[np.searchsorted(tot_light,0.5)]
    return r_half # = r_half/r_scale = r_half/1.0
    
# rs = np.linspace(0,40,100)
# plt.figure()
# plt.plot(rs,sersic(rs,1.5))
# plt.plot(rs,sersic(rs,3.5))

# plt.figure()
# plt.plot(rs,sersic_r(rs,1.5))
# plt.plot(rs,sersic_r(rs,3.5))

# plt.figure()
# plt.plot(rs,sersic_inter(rs,1.5))
# plt.plot(rs,sersic_inter(rs,3.5))

# print sersic_half_scale_ratio(1.5)
# print sersic_half_scale_ratio(3.5)

# plt.show()
redshift=np.linspace(0,2,1000)
mpc_arcsec = interp1d(redshift,cosmo.kpc_proper_per_arcmin(redshift)/1000.0/60)



center_x = 2.5*60*60 #phi #dec
center_y = 87.5*60*60  #theta #ra

post_size_arcmin = 300
post_size = post_size_arcmin*60 #postage stamp size arcsec
post_slct = (cat['ra']>center_y-post_size/2.0 ) & (cat['ra']<center_y+post_size/2.0 )&(cat['dec']>center_x-post_size/2.0 ) & (cat['dec']<center_x+post_size/2.0 )
Mr_obs_slct = cat['magnitude:SDSS_g:observed']<27
slct = post_slct & Mr_obs_slct
print "Done getting slct"
cat = cat[slct]
print "done selecting"
size = cat.shape[0]
slct_blue = cat['magnitude:SDSS_u:observed'].values-cat['magnitude:SDSS_g:observed'].values<1.36
slct_red  = slct_blue==0
redshift = cat['redshift']
dic={}
dic['x(arcsec)'] = cat['dec'].values-center_x
dic['y(arcsec)'] = cat['ra'].values-center_y
dic['gmag']=cat['magnitude:SDSS_g:observed'].values
dic['rmag']=cat['magnitude:SDSS_r:observed'].values
dic['imag']=cat['magnitude:SDSS_i:observed'].values
dic['zmag']=cat['magnitude:SDSS_z:observed'].values
dic['BTg']=bt_func(cat['diskLuminositiesStellar:SDSS_g:observed'].values,cat['spheroidLuminositiesStellar:SDSS_g:observed'].values)
dic['BTr']=bt_func(cat['diskLuminositiesStellar:SDSS_r:observed'].values,cat['spheroidLuminositiesStellar:SDSS_r:observed'].values)
dic['BTi']=bt_func(cat['diskLuminositiesStellar:SDSS_i:observed'].values,cat['spheroidLuminositiesStellar:SDSS_i:observed'].values)
dic['BTz']=bt_func(cat['diskLuminositiesStellar:SDSS_z:observed'].values,cat['spheroidLuminositiesStellar:SDSS_z:observed'].values)
sersic = np.zeros(size)
sersic[slct_blue]=1.5
sersic[slct_red]=3.5
dic['sersic(bulge)']=sersic
half_sphere = np.zeros(size)
half_sphere[slct_blue]=cat['spheroidRadius'][slct_blue].values/sersic_half_scale_ratio(1.5)
half_sphere[slct_red] =cat['spheroidRadius'][slct_red].values/sersic_half_scale_ratio(3.5)
dic['size(Re_bulge)']=half_sphere/mpc_arcsec(redshift)
dic['scale(disklength)']=cat['diskRadius'].values/mpc_arcsec(redshift)
dic['PA']=rand(size)*180.0
ellip=np.zeros(size)
theta = np.arccos(1-2*rand(slct_blue.sum()))
tilt  = (np.cos(theta)+0.1*np.sin(theta))
major = np.maximum(cat['diskRadius'][slct_blue].values,cat['spheroidRadius'][slct_blue].values)
minor = np.maximum(cat['spheroidRadius'][slct_blue].values,tilt*cat['diskRadius'][slct_blue].values)
swap_slct = minor>major
tmp = minor[swap_slct]
minor[swap_slct]=major[swap_slct]
major[swap_slct]=tmp
print major
print minor
print tilt
print minor/major
ellip[slct_blue]=minor/major
ellip[slct_red]=rand(slct_red.sum())*(0.95-0.7)+0.7
dic['ellip']=ellip


def save_dic(file_name,dic):
    print "Saving...."
    keys = ['x(arcsec)','y(arcsec)','gmag','rmag','imag','zmag','sersic(bulge)','BTg','BTr','BTi','BTz','size(Re_bulge)','scale(disklength)','PA','ellip']
    header = ""
    data = []
    spacer = "" #tabs infront of each column name except in the first one
    for key in keys:
        header = header+spacer+key
        spacer = "\t"
        data.append(dic[key])
    np.savetxt(file_name,np.transpose(data),header=header,delimiter='\t',fmt='%.5f',comments="")


save_dic("output/gal_sim/test_full_%ix%i.txt"%(post_size_arcmin,post_size_arcmin),dic)
