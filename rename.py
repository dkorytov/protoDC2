#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import h5py
import dtk
import shutil
import sys
import datetime
import time
import myutil
#from pathlib import Path
from astropy.cosmology import WMAP7 as cosmo
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d 


def swap_slct(slct,val1,val2):
    """swaps the selected rows in val1 with val2"""
    tmp = val1[slct]
    val1[slct] = val2[slct]
    val2[slct] = tmp

def rot_slct(slct,val1,val2,val3):
    """moves the selected rows in val1 into val2, val2 into val3, and val3
    into val1
    """
    tmp = val1[slct]
    val1[slct] = val3[slct]
    val3[slct] = val2[slct]
    val2[slct] = tmp

start_time = time.time()
cosmoAQ = FlatLambdaCDM(H0=71.0,
                        Om0=0.2648)
                       

param = dtk.Param(sys.argv[1])
#output = param.get_string("output")
use_mr_cut = param.get_bool('use_mr_cut')
mr_cut = param.get_int('mr_cut')
gltcs_file_list = param.get_string_list('gltcs_file_list')
lc_rot_info_loc = param.get_string('lc_rot_info_loc')
box = param.get_bool('box')

stepz = dtk.StepZ(200,0,500)
zs = np.linspace(0,3.5,1000)
z_to_dl = interp1d(zs,cosmo.luminosity_distance(zs))

#z_to_cd = interp1d(zs,cosmoAQ.comoving_distance(zs))
output_final = myutil.get_output_loc(param)
gltcs_file = gltcs_file_list[0]
gltcs_param_group = h5py.File(gltcs_file,'r')["Parameters"]
output_mod = output_final.replace('.hdf5','_mod.hdf5')
if(output_final == output_mod):
    print "Can't have the source and modified file be the same!!"
    raise;
print output_mod
hfile = h5py.File(output_mod,'r+')

#################################
### Zeroing out negative masses & their luminosities
#################################

print "Zeroing out..."
comps = ['disk','spheroid']
lum_keys = ['Luminosities','Luminosity','SED','MassStellar','Abundances']
for comp in comps:
    stellar_mass = hfile[comp+'MassStellar'].value
    slct_neg = stellar_mass <= 0;
    for key in hfile.keys():
        if(comp in key and any( lk in key for lk in lum_keys)):
            print '\tzeroing out key',key
            hfile[key][slct_neg] = 0;
print "done zeroing out"



#################################
### Calculating mags and actual observed magnitude 
#################################
print "adding missing mags"

redshifts = hfile['redshift'].value
#cd = z_to_cd(redshifts)
dl = z_to_dl(redshifts)
adjust_mag = -2.5*np.log10(1.0+redshifts)+5*np.log10(dl)+25.0
#adjust_lum = (1.0+redshifts)/dl**2*1e-10
keys = hfile.keys()
#Make the total = disk + spheroid 
sum_vars = ['Luminosity','Luminosities','MassStellar','StarFormationRate','Abundances']
for key in keys:
    if ("disk" in key) and any(sum_var in key for sum_var in sum_vars):
        key_base = key.replace("disk","",1)
        print 'making:  total'+key_base
        hfile["total"+key_base]=hfile["disk"+key_base].value + hfile["spheroid"+key_base].value

keys = hfile.keys()
#Calculate observed magnitude
for key in keys:
    if("totalLuminositiesStellar" in key and  ":observed" in key and ("SDSS" in key or "LSST" in key)):
        key_base = key.replace("totalLuminositiesStellar","",1)
        new_key  = 'magnitude'+key_base
        print "making: "+new_key
        hfile[new_key]=adjust_mag -2.5*np.log10(hfile[key].value)

        # print "making: "+
        # if(new_key not in keys):
            # hfile[new_key]=adjust_mag -2.5*np.log10(hfile[key].value) # 

#Calculate rest magnitude
for key in keys:
    if("totalLuminositiesStellar" in key and  ":rest" in key and ("SDSS" in key or "LSST" in key)):
        key_base = key.replace("totalLuminositiesStellar","",1)
        new_key ="magnitude"+key_base
        print "making: "+new_key
        hfile[new_key]=-2.5*np.log10(hfile[key].value)

lum_r_disk     = hfile['diskLuminositiesStellar:SDSS_r:rest'].value
lum_r_spheroid = hfile['spheroidLuminositiesStellar:SDSS_r:rest'].value
lum_r_tot = lum_r_disk+lum_r_spheroid
# ellip = hfile['diskEllipticity'].value*lum_r_disk/lum_r_tot + hfile['spheroidEllipticity'].value*lum_r_spheroid/lum_r_tot
# hfile['totalEllipticity']=ellip
# pos_angle = hfile['positionAngle'].value
# print "pos_angle",np.min(pos_angle),np.max(pos_angle)
# hfile['totalEllipticity1'] = ellip*np.cos(pos_angle*2.0)
# hfile['totalEllipticity2'] = ellip*np.sin(pos_angle*2.0)

print hfile.keys()
mag_r = hfile['magnitude:SDSS_r:rest'].value
sm = hfile['totalMassStellar'].value
print "print min max mag r", np.min(mag_r),np.max(mag_r)
print "print min max stellar mass", np.min(sm), np.max(sm)


# redshifts = result['redshift']
# print "min redshift: ",np.min(redshifts)
# print "max redshift: ",np.max(redshifts)
# dl = z_to_dl(redshifts)
# print "done getting the luminosity distance..."
# adjust = -2.5*np.log10(1+redshifts)+5*np.log10(dl)+25.0
# for key in gal_df.keys():
#     print 'test observed mag', key
#     if("diskLuminositiesStellar" in key and  ":observed" in key and ("SDSS" in key or "LSST" in key)):
#         print '\t pass'
#         key_base = key.replace("diskLuminositiesStellar","")
#         result['magnitude'+key_base]=adjust -2.5*np.log10(result['diskLuminositiesStellar'+key_base]+result['spheroidLuminositiesStellar'+key_base])

    
print "\tadded. "


#################################
###Add units to all fields
#################################
mag_list = ['magnitude']; mag_unit = 'AB magnitude'
arcsec_list= ['Arcsec']; arcsec_unit = 'arcsecond'
rad_list = []; rad_unit ='radians'
deg_list = ['ra','dec','positionAngle','inclination']; deg_unit = 'degrees'
phys_kpc_list = ['Radius']; phys_kpc_unit = 'physical kpc'
phys_mpc_list = []; phys_mpc_unit = 'physical Mpc'
reduced_dist_list =['Reduced','EigenVector'];reduced_dist_unit = 'unitless'
eigen_val_list = ['EigenValue'];eigen_val_unit = 'comoving Mpc/h'
comv_mpc_list = ['x','y','z']; comv_mpc_unit = 'comoving Mpc/h'
vel_list = ['vx','vy','vz','Velocity']; vel_unit = 'km/s'
timeSFR_list =['TimeWeightedIntegratedSFR']; timeSFR_unit = 'Gyr*Msun'
sfr_list =['SFR','blackHoleAccretionRate','StarFormationRate']; sfr_unit = 'Msun/Gyr'
mass_list =['Mass','IntegratedSFR']; mass_unit = 'Msun'
abundance_list =['Abundance'];abundance_unit = 'Msun'
luminosity_list =['Luminosities','Luminosity']; luminosity_unit = 'AB luminosity (4.4659e13 W/Hz)'
unitless_list = ['redshift','shear','magnification','convergence','Ellipticity','Sersic','AxisRatio']; unitless_unit ='unitless'
id_list = ['Index','Tag','placementType','galaxyID']; id_unit = 'id/index'
angular_list = ['angularMomentum'];angular_unit = 'Msun*km/s*Mpc'
bool_list =['nodeIsIsolated'];bool_unit = 'boolean'
spinSpin_list =['spinSpin'];spinSpin_unit ='lambda'
step_list = ['step'];step_unit = 'simluation step'
print "assigning units"
for key in hfile.keys():
    #add an empty discription
    # desc_file_loc = Path(desc_loc+'/'+key)
    # if desc_file_loc.is_file():
    #     with open(desc_file_loc) as desc_file:
    #         hfile[key].attrs['description']=desc_file.read()
    # else:
    #     with open(desc_file_loc) as desc_file:
    #         desc_file.write("None given")
    #     hfile[key].attrs['description']='None given'
    print key
    print any(l in key for l in comv_mpc_list)
    #add magnitude units
    if(any(l in key for l in mag_list)):
        hfile[key].attrs['units']=mag_unit
        print "\t mag"
    #add arcsec units
    elif(any(l in key for l in arcsec_list)):
        hfile[key].attrs['units']=arcsec_unit
        print "\t ",arcsec_unit
    #add rad units
    elif(any(l in key for l in rad_list)):
        hfile[key].attrs['units']=rad_unit
        print "\t ",rad_unit
    #add degree units
    elif(any(l in key for l in deg_list)):
        hfile[key].attrs['units']=deg_unit
        print '\t',deg_unit
    #add kpc units
    elif(any(l in key for l in phys_kpc_list)):
        hfile[key].attrs['units']=phys_kpc_unit
        print "\t ",phys_kpc_unit
    #add mpc units
    elif(any(l in key for l in phys_mpc_list)):
        hfile[key].attrs['units']=phys_mpc_unit
        print "\t ",phys_mpc_unit
    #reduced distances units
    elif(any(l in key for l in reduced_dist_list)):
        hfile[key].attrs['units']=reduced_dist_unit
        print "\t ",reduced_dist_unit
    #eigen val units
    elif(any(l in key for l in eigen_val_list)):
        hfile[key].attrs['units']=eigen_val_unit
        print "\t ",reduced_dist_unit
    #add comoving mpc units
    elif(any(l == key for l in comv_mpc_list)):
        hfile[key].attrs['units']=comv_mpc_unit
        print "\t ",comv_mpc_unit
    #add velocity units
    elif(any(l in key for l in vel_list)):
        hfile[key].attrs['units']=vel_unit
        print "\t ",vel_unit
    #add timesfr
    elif(any(l in key for l in timeSFR_list)):
        hfile[key].attrs['units']=timeSFR_unit
        print "\t ",timeSFR_unit
    #add sfr
    elif(any(l in key for l in sfr_list)):
        hfile[key].attrs['units']=sfr_unit
        print "\t ",sfr_unit
    #add mass
    elif(any(l in key for l in mass_list)):
        hfile[key].attrs['units']=mass_unit
        print "\t ",mass_unit
    #add abundance
    elif(any(l in key for l in abundance_list)):
        hfile[key].attrs['units']=abundance_unit
        print "\t ",abundance_unit
    #add luminosity units
    elif(any(l in key for l in luminosity_list)):
        hfile[key].attrs['units']=luminosity_unit
        print "\t ",luminosity_unit
    #add unit less
    elif(any(l in key for l in unitless_list)):
        hfile[key].attrs['units']=unitless_unit
        print "\t ",unitless_unit
    #add mass units
    elif(any(l in key for l in id_list)):
        hfile[key].attrs['units']=id_unit
        print "\t ",id_unit
    #angular momentum 
    elif(any(l in key for l in angular_list)):
        hfile[key].attrs['units']=angular_unit
        print "\t ",angular_unit
    #boolean
    elif(any(l in key for l in bool_list)):
        hfile[key].attrs['units']=bool_unit
        print "\t", bool_unit
    #spinSpin
    elif(any(l in key for l in spinSpin_list)):
        hfile[key].attrs['units']=spinSpin_unit
    #step
    elif(any(l in key for l in step_list)):
        hfile[key].attrs['units']=step_unit
    #Everything should have a unit!
    else:
        print "column", key, "was not assigned a unit :("
        raise;
print "done assigning"

#################################
print "renaming..."
#rename nodeIsIsolated
hfile.move('nodeIsIsolated','isCentral')

#remove parentIndex
#del hfile['parentIndex']

#rename nodeIndex
hfile.move('nodeIndex','infallIndex')
#shift ra/dec & swap lensed & true
if not box:
    hfile.move('ra','ra_true')
    hfile.move('dec','dec_true')
    hfile.move('ra_lensed','ra')
    hfile.move('dec_lensed','dec')
#rename position type
# position_type = hfile['placementType'].value
# slct_central = position_type == 5
# slct_core    = position_type == 2
# slct_infall  = position_type == 0
# slct_rnd     = position_type == 1
# slct_nfw     = position_type == 5
# slct_blank   = position_type == 3
# position_type[slct_central]=0
# position_type[slct_core]=1
# position_type[slct_infall]=2
# position_type[slct_rnd]=3
# position_type[slct_nfw]=4
# position_type[slct_blank]=5
# #TODO check position types

#del hfile['placementType']
# hfile['placementType'][...]=position_type
#TODO check that attrs stay the same


#rename redshifts to have redshift & redshift_hubble
# del hfile['redshift']
# hfile.move('redshiftObserver','redshift')
# hfile.move('redshiftCosmological','redshiftHubble')

#change to degrees
if not box:
    hfile['ra_true'][...] = hfile['ra_true'].value/3600.0-87.5
    hfile['dec_true'][...] = hfile['dec_true'].value/3600.0-2.5
    hfile['ra'][...] = hfile['ra'].value/3600.0
    hfile['dec'][...] = hfile['dec'].value/3600.0
# hfile['positionAngle'][...]=hfile['positionAngle'].value*180.0/np.pi
# hfile.move('diskRadius','diskScaleRadius')
# hfile.move('spheroidRadius','spheroidScaleRadius')
# hfile.move('diskRadiusArcsec','diskScaleRadiusArcsec')
# hfile.move('spheroidRadiusArcsec','spheroidScaleRadiusArcsec')

print "done renaming"


#################################
### Rotating halo shapes
if not box:
    rot_hfile = h5py.File(lc_rot_info_loc,'r')
    gal_id = rot_hfile['galaxyID'].value
    rotations = rot_hfile['rotations'].value
    srt = np.argsort(gal_id)
    indx = dtk.search_sorted(gal_id,hfile['galaxyID'].value,sorter=srt)
    slct_neg = indx == -1
    if(np.sum(slct_neg) > 0):
        #this should never happen. Every galaxy should have a rotation val
        #unless the cuts changed
        print "not matched: ", np.sum(slct_neg),slct_neg.shape
        raise;
    rots = rotations[indx]

    #rots == 0 unchanged
    #1  x <-> y 
    #2  x <-> z
    #3  y <-> z
    #4  x->z y->x z->y
    #5  x->y y->z z->x
    hfile['lightConeRotation']=rots
    hfile['lightConeRotation'].attrs['units']="id/index"
    slct = rots == 1
    for vect_type in ("","Reduced"):
        for i in range(1,4):
            print 'hostHaloEigenVector%s%dX'%(vect_type,i)
            swap_slct(slct,
                      hfile['hostHaloEigenVector%s%dX'%(vect_type,i)],
                      hfile['hostHaloEigenVector%s%dY'%(vect_type,i)])

    slct = rots == 2
    for vect_type in ("","Reduced"):
        for i in range(1,4):
            swap_slct(slct,
                      hfile['hostHaloEigenVector%s%dX'%(vect_type,i)],
                      hfile['hostHaloEigenVector%s%dZ'%(vect_type,i)])
    slct = rots == 3
    for vect_type in ("","Reduced"):
        for i in range(1,4):
            swap_slct(slct,
                      hfile['hostHaloEigenVector%s%dY'%(vect_type,i)],
                      hfile['hostHaloEigenVector%s%dZ'%(vect_type,i)])

    slct = rots == 4
    for vect_type in ("","Reduced"):
        for i in range(1,4):
            rot_slct(slct,
                     hfile['hostHaloEigenVector%s%dX'%(vect_type,i)],
                     hfile['hostHaloEigenVector%s%dY'%(vect_type,i)],
                     hfile['hostHaloEigenVector%s%dZ'%(vect_type,i)])

    slct = rots == 5
    for vect_type in ("","Reduced"):
        for i in range(1,4):
            rot_slct(slct,
                     hfile['hostHaloEigenVector%s%dZ'%(vect_type,i)],
                     hfile['hostHaloEigenVector%s%dY'%(vect_type,i)],
                     hfile['hostHaloEigenVector%s%dX'%(vect_type,i)])


################################
###Create groups
#galaxyProp
print "moving to groups"
morphology_list = ["Radius", "Major", "Minor", "Ellipticity", "inclination",'positionAngle','AxisRatio','Sersic']
for key in hfile.keys():
    print '\t\t',key
    if("SED" in key):
        hfile.move(key,'galaxyProperties/SEDs/'+key)
    elif("SDSS" in key):
        hfile.move(key,'galaxyProperties/SDSS_filters/'+key)
    elif("LSST" in key):
        hfile.move(key,'galaxyProperties/LSST_filters/'+key)
    elif("Line" in key):
        hfile.move(key,'galaxyProperties/emissionLines/'+key)
    elif("Luminosities" in key or "Luminosity" in key):
        hfile.move(key,'galaxyProperties/otherLuminosities/'+key)
    elif(any(l in key for l in morphology_list)):
        hfile.move(key,'galaxyProperties/morphology/'+key)
    else:
        hfile.move(key,'galaxyProperties/'+key)






################################
### Add meta data
print "adding metadata"

hgroup = hfile.create_group('metaData')
hgroup['ellipticityDefinition']='''(a-b)/(a+b), a = major axis, b = minor axis
ellipticity1 = ellipticity * cos(2 alpha)
ellipticity2 = ellipticity * sin(2 alpha)
where alpha = position angle
'''
hgroup['shearDefinition']='a11 = 1-k-gamma1, a22 =1-k+gamma1 , a21=a12=gamma2'

hgroup['shapeModel']="""Spheroid Componet: 
The profile is assumed to be Hernquist profile which has a Sersic
index = 4.  The axis ratio is a uniform distribution from 0.7 to 0.95

Disk Componet:

The profile is assumed to Sersic index =1. The disk thickness is
assumed to be 10% of the disk radius. For the 'b' major axis, 'a' minor
axis, 'R' scale radius and 'theta' inclination:
a = R
b = R*sqrt(cos(theta)**2 + (0.1*sin(theta))**2)
Note: for small theta, a is larger than b. In that case, they are swapped.

The positionAngle is uniform random angle from 0 -> pi. 
"""
hfile['/metaData/skyArea']=25.0

hgroup = hfile.create_group('metaData/simulationParameters')
hgroup['H_0']=71.0
hgroup['H_0'].attrs['units']='km/s/Mpc'
hgroup['Omega_DE']=0.7352
hgroup['Omega_matter']=0.2648
hgroup['Omega_b']=0.0448
hgroup['Omega_Nu']=0.0
hgroup['w_de']=-1.0
hgroup['sigma_8']=0.8
hgroup['N_s']=0.963
hgroup['boxSize']=360.563
hgroup['boxSize'].attrs['units']='comoving Mpc'
hgroup['NP']=1024
hgroup['particleMass']=1.3709e9/0.71
hgroup['particleMass'].attrs['units']='Msun'
hgroup['haloMassDefinition']='FoF b=0.168'
#group = hfile.create_group('metaData/GalacticusParameters')
hfile['/metaData/catalogCreationDate']=datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
cut_str = ""
if(use_mr_cut):
    cut_str += 'magnitude:SDSS_r:rest < %.2f'%mr_cut;
else:
    cut_str += 'no cuts'
hfile['/metaData/cuts']=cut_str

hfile['/metaData/versionChangeNotes'] = """version notes:
1) Correct inclinations bug (degrees =/= radians)
2) Changed the model for the inclination effect on the minor axis
3) Added shape model description in the metaData
4) Changed ellipticity definition from e=1-q to e = (1-q)/(1+q)
5) Added ellipticity_1 & 2
6) Added ellipticity definition to metaData 
7) Added axis ratios data columns
8) Changed Ra/Dec, inclination, and position angle to be in degrees. 
9) Any disk/spheroid component with zero or less stellar mass, has all related 
luminosites set to zero.
10) Shear bug fix
"""
hfile['/metaData/version']='2.1.2'
hfile['/metaData/versionMajor']=2
hfile['/metaData/versionMinor']=1
hfile['/metaData/versionMinorMinor']=2


hfile.copy(gltcs_param_group,'/metaData/GalacticusParameters')

print "done", time.time()-start_time
#describe position type
