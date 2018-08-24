#!/usr/bin/env python2.7

import numpy as np
import h5py 
import datetime
from astropy.cosmology import WMAP7 as cosmo
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from numpy.random import rand
from numpy.random import seed
from pecZ import *
import dtk
import sys 

def multi_slct_assign(dest, data_in,slct1,slct2):
    slct = sclt1
    slct[slct1]=slct2
    dest[sclt1]=data_in

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

param = dtk.Param(sys.argv[1])
output = param.get_string('output')
sod_loc = param.get_string("sod_loc")
halo_shape_loc = param.get_string("halo_shape_loc")
halo_shape_red_loc = param.get_string("halo_shape_red_loc")
rnd_seed = param.get_int("rnd_seed")
seed(rnd_seed)
hfile = h5py.File(output,'a')
hfile.attrs['H_0']=71.0
hfile.attrs['Omega_DE']=0.7352
hfile.attrs['Omega_matter']=0.2648
hfile.attrs['Omega_b']=0.0448
hfile.attrs['Omega_Nu']=0.0
hfile.attrs['w_de']=-1.0
hfile.attrs['sigma_8']=0.8
hfile.attrs['N_s']=0.963
hfile.attrs['box_size']=360.5636
hfile.attrs['NP']=1024
hfile.attrs['particle_mass']=1.3709e9/0.71
hfile.attrs['catalog_creation_date']=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

if('theta' in hfile.keys()):
    hfile.move('theta','ra')
if('phi' in hfile.keys()):
    hfile.move('phi','dec')
if('theta_lensed' in hfile.keys()):
    hfile.move('theta_lensed','ra_lensed')
if('phi_lensed' in hfile.keys()):
    hfile.move('phi_lensed','dec_lensed')

if('k0' in hfile.keys()):
    hfile.move('k0','convergence')
if('m0' in hfile.keys()):
    hfile.move('m0','magnification')
if('galaxyID' not in hfile.keys()):
    length = hfile['nodeIndex'].len()
    print length
    hfile['galaxyID']=np.arange(0,length,dtype='i8')

try:
    del hfile['redshiftObserver']
    del hfile['redshiftCosmological']
    del hfile['peculiarVelocity']
except KeyError:
    print("redshift values haven't been calculated yet, all is fine")

if 'x' in hfile.keys():
    redshift = hfile['redshift'].value
    x = hfile['x'].value
    y = hfile['y'].value
    z = hfile['z'].value
    vx = hfile['vx'].value
    vy = hfile['vy'].value
    vz = hfile['vz'].value
    z_hubb = hfile['redshift'].value
    z_pec,z_tot,v_pec,v_pec_a,r1,r2,r_dist = pecZ(x,y,z,vx,vy,vz,redshift)
    hfile['redshiftObserver'] = z_tot
    hfile['redshiftCosmological'] = redshift
    hfile['peculiarVelocity'] = v_pec



#Adding size as half light radius
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


redshift=np.linspace(0,2,1000)
#mpc_arcsec = interp1d(redshift,cosmo.kpc_proper_per_arcmin(redshift)/1000.0/60)
arcsec_per_mpc = interp1d(redshift,cosmo.arcsec_per_kpc_proper(redshift).value*1000.0)
sersic_half_scale_ratio_10 = sersic_half_scale_ratio(1.0)
sersic_half_scale_ratio_15 = sersic_half_scale_ratio(1.5)
sersic_half_scale_ratio_35 = sersic_half_scale_ratio(3.5)

hernquist_half_scale_ratio = 1.8153
size = hfile['redshift'].size

# slct_blue = hfile['magnitude:SDSS_u:observed'].values-hfile['magnitude:SDSS_g:observed'].values<1.36
# slct_red  = slct_blue==0
# size = slct_blue.size
# data = np.ndarray(slct_blue.size,dtype=float)

# data[slct_blue]=hfile['spheriodRadius'].values*sersic_half_scale_ratio_15
# data[slct_red] =hfile['spheriodRadius'].values*sersic_half_scale_ratio_35
def delhfile(hfile,var):
    if var in hfile.keys():
        del hfile[var]
try:
    print "deleting old versions of variables"
    delhfile( hfile,'spheroidRadiusHalfLight')
    delhfile( hfile,'diskRadiusHalfLight') 
    delhfile( hfile,'spheroidRadiusArcsec')
    delhfile( hfile,'spheroidRadiusHalfLightArcsec')
    delhfile( hfile,'diskRadiusArcsec')
    delhfile( hfile,'diskRadiusHalfLightArcsec')
    delhfile( hfile,'spheroidHalfLightRadius')
    delhfile( hfile,'diskHalfLightRadius') 
    delhfile( hfile,'spheroidScaleRadiusArcsec')
    delhfile( hfile,'spheroidHalfLightRadiusArcsec')
    delhfile( hfile,'diskScaleRadiusArcsec')
    delhfile( hfile,'diskHalfLightRadiusArcsec')
    delhfile( hfile,'diskEllipticity')
    delhfile( hfile,'spheroidEllipticity')
    delhfile( hfile,'alignment')
    delhfile( hfile,'diskMajorAxisArcsec')
    delhfile( hfile,'diskMinorAxisArcsec')
    delhfile( hfile,'spheroidMinorAxisArcsec')
    delhfile( hfile,'spheroidMajorAxisArcsec')
    delhfile( hfile,'spheroidSersicIndex')
    delhfile( hfile,'diskSersicIndex')
    delhfile( hfile,'spheroidSersicIndex')
    delhfile( hfile,'spheroidAxisRatio')
    delhfile( hfile,'diskAxisRatio')
    delhfile( hfile,'positionAngle')
    delhfile( hfile,'spheroidEllipticity1')
    delhfile( hfile,'spheroidEllipticity2')
    delhfile( hfile,'diskEllipticity1')
    delhfile( hfile,'diskEllipticity2')
    delhfile( hfile,'hostHaloEigenValue1')
    delhfile( hfile,'hostHaloEigenValue2')
    delhfile( hfile,'hostHaloEigenValue3')
    delhfile( hfile,'hostHaloEigenVector1X')
    delhfile( hfile,'hostHaloEigenVector1Y')
    delhfile( hfile,'hostHaloEigenVector1Z')
    delhfile( hfile,'hostHaloEigenVector2X')
    delhfile( hfile,'hostHaloEigenVector2Y')
    delhfile( hfile,'hostHaloEigenVector2Z')
    delhfile( hfile,'hostHaloEigenVector3X')
    delhfile( hfile,'hostHaloEigenVector3Y')
    delhfile( hfile,'hostHaloEigenVector3Z')
    delhfile( hfile,'hostHaloEigenValueReduced1')
    delhfile( hfile,'hostHaloEigenValueReduced2')
    delhfile( hfile,'hostHaloEigenValueReduced3')
    delhfile( hfile,'hostHaloEigenVectorReduced1X')
    delhfile( hfile,'hostHaloEigenVectorReduced1Y')
    delhfile( hfile,'hostHaloEigenVectorReduced1Z')
    delhfile( hfile,'hostHaloEigenVectorReduced2X')
    delhfile( hfile,'hostHaloEigenVectorReduced2Y')
    delhfile( hfile,'hostHaloEigenVectorReduced2Z')
    delhfile( hfile,'hostHaloEigenVectorReduced3X')
    delhfile( hfile,'hostHaloEigenVectorReduced3Y')
    delhfile( hfile,'hostHaloEigenVectorReduced3Z')
    delhfile( hfile,'hostHaloSODTag')
    delhfile( hfile,'hostHaloSODMass')
except KeyError:
    print("can't delete them all...could be fine, this is just a overwrite work around2")
hfile['positionAngle']=np.array(np.pi*rand(size))
hfile['spheroidSersicIndex']=4*np.ones(hfile['redshift'].value.size)
hfile['diskSersicIndex']=1.0*np.ones(hfile['redshift'].value.size)
data = hfile['spheroidRadius'].value*hernquist_half_scale_ratio
hfile['spheroidHalfLightRadius']=data
hfile['diskHalfLightRadius'] = hfile['diskRadius'].value*sersic_half_scale_ratio_10
#gal_mpc_to_arcsec = 1.0/mpc_arcsec(hfile['redshift'].value)
gal_mpc_to_arcsec = arcsec_per_mpc(hfile['redshift'].value)
hfile['spheroidRadiusArcsec']=hfile['spheroidRadius'].value*gal_mpc_to_arcsec
hfile['spheroidHalfLightRadiusArcsec']=hfile['spheroidHalfLightRadius'].value*gal_mpc_to_arcsec
hfile['diskRadiusArcsec']=hfile['diskRadius'].value*gal_mpc_to_arcsec
hfile['diskHalfLightRadiusArcsec']=hfile['diskRadius'].value*gal_mpc_to_arcsec
q = 0.7+rand(size)*(0.95-0.7)
q = dtk.clipped_gaussian(0.8, 0.2, size, max_val = 1.0, min_val=0.0)
hfile['spheroidAxisRatio'] = q
hfile['spheroidMajorAxisArcsec'] = hfile['spheroidHalfLightRadiusArcsec'].value
hfile['spheroidMinorAxisArcsec'] = hfile['spheroidHalfLightRadiusArcsec'].value*q
hfile['spheroidEllipticity'] = (1.0-q)/(1.0+q)
hfile['spheroidEllipticity1'] = np.cos(2.0*hfile['positionAngle'].value)*hfile['spheroidEllipticity'].value
hfile['spheroidEllipticity2'] = np.sin(2.0*hfile['positionAngle'].value)*hfile['spheroidEllipticity'].value
#theta = np.arccos(1-2*rand(size))
print "Thetas"
print hfile['inclination'].value
theta = hfile['inclination'].value*np.pi/180.0
dist,lim = dtk.make_distribution(-theta)
resamp = dtk.resample_distribution(dist,gal_zoo_dist,lim,[0.0,1.0])

tilt = resamp(-theta)
print "title: ",np.min(tilt),np.max(tilt)
#tilt = np.sqrt(np.cos(theta)**2+0.01*np.sin(theta)**2)
major = hfile['diskHalfLightRadiusArcsec']
minor = hfile['diskHalfLightRadiusArcsec']*tilt
#make sure that major is infact larger
# swap_slct = minor>major
# tmp = minor[swap_slct]
# minor[swap_slct]=major[swap_slct]
# major[swap_slct]=tmp
print hfile.keys()
hfile['diskMinorAxisArcsec'] = minor
hfile['diskMajorAxisArcsec'] = major
q = minor/major
hfile['diskAxisRatio']=q
hfile['diskEllipticity']=(1.0-q)/(1.0+q)
hfile['diskEllipticity1'] = np.cos(2.0*hfile['positionAngle'].value)*hfile['diskEllipticity'].value
hfile['diskEllipticity2'] = np.sin(2.0*hfile['positionAngle'].value)*hfile['diskEllipticity'].value


#################################
### Getting additional host halo info
#################################
steps = hfile['step'].value
htag = hfile['hostHaloTag'].value
htag_real = dtk.frag_to_real(htag)
so_mass = np.zeros_like(steps,dtype='float')
halo_eg1 = np.zeros_like(steps,dtype='float')
halo_eg2 = np.zeros_like(steps,dtype='float')
halo_eg3 = np.zeros_like(steps,dtype='float')
halo_eg1_x = np.zeros_like(steps,dtype='float')
halo_eg1_y = np.zeros_like(steps,dtype='float')
halo_eg1_z = np.zeros_like(steps,dtype='float')
halo_eg2_x = np.zeros_like(steps,dtype='float')
halo_eg2_y = np.zeros_like(steps,dtype='float')
halo_eg2_z = np.zeros_like(steps,dtype='float')
halo_eg3_x = np.zeros_like(steps,dtype='float')
halo_eg3_y = np.zeros_like(steps,dtype='float')
halo_eg3_z = np.zeros_like(steps,dtype='float')

#reduced eigen shapes
halo_egr1 = np.zeros_like(steps,dtype='float')
halo_egr2 = np.zeros_like(steps,dtype='float')
halo_egr3 = np.zeros_like(steps,dtype='float')
halo_egr1_x = np.zeros_like(steps,dtype='float')
halo_egr1_y = np.zeros_like(steps,dtype='float')
halo_egr1_z = np.zeros_like(steps,dtype='float')
halo_egr2_x = np.zeros_like(steps,dtype='float')
halo_egr2_y = np.zeros_like(steps,dtype='float')
halo_egr2_z = np.zeros_like(steps,dtype='float')
halo_egr3_x = np.zeros_like(steps,dtype='float')
halo_egr3_y = np.zeros_like(steps,dtype='float')
halo_egr3_z = np.zeros_like(steps,dtype='float')

print steps
steps_unique = np.unique(steps)
print "loading in halo shapes and extra halo info"

for step in steps_unique:
    print "working on step: ", step
    slct_step = steps == step
    print step,sod_loc
    sod_step_loc = sod_loc.replace('${step}',str(step))
    print sod_step_loc
    halo_shape_step_loc = halo_shape_loc.replace('${step}',str(step))
    halo_shape_red_step_loc = halo_shape_red_loc.replace('${step}',str(step))
    sod_cat_htag = dtk.gio_read(sod_step_loc,'fof_halo_tag')
    sod_cat_mass = dtk.gio_read(sod_step_loc,'sod_halo_mass')
    srt = np.argsort(sod_cat_htag)
    indx = dtk.search_sorted(sod_cat_htag,htag_real[slct_step],sorter=srt)
    slct_indx = indx != -1
    slct = slct_step
    slct[slct_step]=slct_indx
    #Assign so masses to the ones that have it
    so_mass[slct] = sod_cat_mass[indx[slct_indx]]
    #Assign a dumby value to those that don't have it
    #slct[slct_step]=(slct_indx==0)
    #so_mass[slct] = -1.0
    eg_cat_htag = dtk.gio_read(halo_shape_step_loc,'halo_id')
    eg_cat_eg1 = dtk.gio_read(halo_shape_step_loc,'eval1')
    eg_cat_eg2 = dtk.gio_read(halo_shape_step_loc,'eval2')
    eg_cat_eg3 = dtk.gio_read(halo_shape_step_loc,'eval3')
    eg_cat_eg1_x = dtk.gio_read(halo_shape_step_loc,'evec1x')
    eg_cat_eg1_y = dtk.gio_read(halo_shape_step_loc,'evec1y')
    eg_cat_eg1_z = dtk.gio_read(halo_shape_step_loc,'evec1z')
    eg_cat_eg2_x = dtk.gio_read(halo_shape_step_loc,'evec2x')
    eg_cat_eg2_y = dtk.gio_read(halo_shape_step_loc,'evec2y')
    eg_cat_eg2_z = dtk.gio_read(halo_shape_step_loc,'evec2z')
    eg_cat_eg3_x = dtk.gio_read(halo_shape_step_loc,'evec3x')
    eg_cat_eg3_y = dtk.gio_read(halo_shape_step_loc,'evec3y')
    eg_cat_eg3_z = dtk.gio_read(halo_shape_step_loc,'evec3z')
    srt = np.argsort(eg_cat_htag)
    indx = dtk.search_sorted(eg_cat_htag,htag_real[slct_step],sorter=srt)
    slct_indx = indx != -1
    slct = slct_step
    slct[slct_step]=slct_indx
    #assign value to those who have it
    halo_eg1[slct] = eg_cat_eg1[indx[slct_indx]]
    halo_eg2[slct] = eg_cat_eg2[indx[slct_indx]]
    halo_eg3[slct] = eg_cat_eg3[indx[slct_indx]]
    halo_eg1_x[slct] = eg_cat_eg1_x[indx[slct_indx]]
    halo_eg1_y[slct] = eg_cat_eg1_y[indx[slct_indx]]
    halo_eg1_z[slct] = eg_cat_eg1_z[indx[slct_indx]]
    halo_eg2_x[slct] = eg_cat_eg2_x[indx[slct_indx]]
    halo_eg2_y[slct] = eg_cat_eg2_y[indx[slct_indx]]
    halo_eg2_z[slct] = eg_cat_eg2_z[indx[slct_indx]]
    halo_eg3_x[slct] = eg_cat_eg3_x[indx[slct_indx]]
    halo_eg3_y[slct] = eg_cat_eg3_y[indx[slct_indx]]
    halo_eg3_z[slct] = eg_cat_eg3_z[indx[slct_indx]]
    
    #Now doing the reduced values
    eg_cat_htag = dtk.gio_read(halo_shape_red_step_loc,'halo_id')
    eg_cat_eg1 = dtk.gio_read(halo_shape_red_step_loc,'eval1')
    eg_cat_eg2 = dtk.gio_read(halo_shape_red_step_loc,'eval2')
    eg_cat_eg3 = dtk.gio_read(halo_shape_red_step_loc,'eval3')
    eg_cat_eg1_x = dtk.gio_read(halo_shape_red_step_loc,'evec1x')
    eg_cat_eg1_y = dtk.gio_read(halo_shape_red_step_loc,'evec1y')
    eg_cat_eg1_z = dtk.gio_read(halo_shape_red_step_loc,'evec1z')
    eg_cat_eg2_x = dtk.gio_read(halo_shape_red_step_loc,'evec2x')
    eg_cat_eg2_y = dtk.gio_read(halo_shape_red_step_loc,'evec2y')
    eg_cat_eg2_z = dtk.gio_read(halo_shape_red_step_loc,'evec2z')
    eg_cat_eg3_x = dtk.gio_read(halo_shape_red_step_loc,'evec3x')
    eg_cat_eg3_y = dtk.gio_read(halo_shape_red_step_loc,'evec3y')
    eg_cat_eg3_z = dtk.gio_read(halo_shape_red_step_loc,'evec3z')
    srt = np.argsort(eg_cat_htag)
    indx = dtk.search_sorted(eg_cat_htag,htag_real[slct_step],sorter=srt)
    slct_indx = indx != -1
    slct = slct_step
    slct[slct_step]=slct_indx
    halo_egr1[slct] = eg_cat_eg1[indx[slct_indx]]
    halo_egr2[slct] = eg_cat_eg2[indx[slct_indx]]
    halo_egr3[slct] = eg_cat_eg3[indx[slct_indx]]
    halo_egr1_x[slct] = eg_cat_eg1_x[indx[slct_indx]]
    halo_egr1_y[slct] = eg_cat_eg1_y[indx[slct_indx]]
    halo_egr1_z[slct] = eg_cat_eg1_z[indx[slct_indx]]
    halo_egr2_x[slct] = eg_cat_eg2_x[indx[slct_indx]]
    halo_egr2_y[slct] = eg_cat_eg2_y[indx[slct_indx]]
    halo_egr2_z[slct] = eg_cat_eg2_z[indx[slct_indx]]
    halo_egr3_x[slct] = eg_cat_eg3_x[indx[slct_indx]]
    halo_egr3_y[slct] = eg_cat_eg3_y[indx[slct_indx]]
    halo_egr3_z[slct] = eg_cat_eg3_z[indx[slct_indx]]


hfile['hostHaloSODTag'] = htag_real
hfile['hostHaloSODMass'] = so_mass

hfile['hostHaloEigenValue1']=halo_eg1
hfile['hostHaloEigenValue2']=halo_eg2
hfile['hostHaloEigenValue3']=halo_eg3
hfile['hostHaloEigenVector1X']=halo_eg1_x
hfile['hostHaloEigenVector1Y']=halo_eg1_y
hfile['hostHaloEigenVector1Z']=halo_eg1_z
hfile['hostHaloEigenVector2X']=halo_eg2_x
hfile['hostHaloEigenVector2Y']=halo_eg2_y
hfile['hostHaloEigenVector2Z']=halo_eg2_z
hfile['hostHaloEigenVector3X']=halo_eg3_x
hfile['hostHaloEigenVector3Y']=halo_eg3_y
hfile['hostHaloEigenVector3Z']=halo_eg3_z

hfile['hostHaloEigenValueReduced1']=halo_eg1
hfile['hostHaloEigenValueReduced2']=halo_eg2
hfile['hostHaloEigenValueReduced3']=halo_eg3
hfile['hostHaloEigenVectorReduced1X']=halo_eg1_x
hfile['hostHaloEigenVectorReduced1Y']=halo_eg1_y
hfile['hostHaloEigenVectorReduced1Z']=halo_eg1_z
hfile['hostHaloEigenVectorReduced2X']=halo_eg2_x
hfile['hostHaloEigenVectorReduced2Y']=halo_eg2_y
hfile['hostHaloEigenVectorReduced2Z']=halo_eg2_z
hfile['hostHaloEigenVectorReduced3X']=halo_eg3_x
hfile['hostHaloEigenVectorReduced3Y']=halo_eg3_y
hfile['hostHaloEigenVectorReduced3Z']=halo_eg3_z


# 'correcting spelling mistakes..'
for key in hfile.keys():
    if ('magnitdue'in key):
        key2 = key.replace('magnitdue','magnitude')
        print key,'->',key2
        hfile.move(key,key2)

hfile.close()
