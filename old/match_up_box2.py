#!/usr/bin/env python2.7


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from astropy.cosmology import WMAP7 as cosmo
import sys
import time
from scipy.interpolate import interp1d 
import dtk


param = dtk.Param(sys.argv[1])
lc_loc = param.get_string("lc_loc")
gltcs_loc = param.get_string("gltcs_loc")
shear_loc = param.get_string("shear_loc")
gio_loc   = param.get_string("gio_loc")
use_shear = param.get_bool("use_shear")
gltcs_file_num = param.get_int("gltcs_file_num")
gltcs_file_step = param.get_int("gltcs_file_step")
gltcs_steps     = param.get_int_list("gltcs_steps")
gltcs_internal  = param.get_int_list("gltcs_internal")
gltcs_str_z     = param.get_string_list("gltcs_str_z")
use_mr_cut      = param.get_bool("use_mr_cut")
mr_cut          = param.get_float("mr_cut")
tmp_to_disk     = param.get_bool("tmp_to_disk")
tmp_file_loc    = param.get_string("tmp_file_loc")
output          = param.get_string("output")
output_box = param.get_string("output_box")
stepz = dtk.StepZ(200,0,500)
stepz = dtk.StepZ(200,0,500)
zs = np.linspace(0,1.5,1000)
z_to_dl = interp1d(zs,cosmo.luminosity_distance(zs))
total_start = time.time()

def catalog_to_hdf5(cat,output):
    hfile = h5py.File(output,'w')
    steps = cat.get_steps()
    var_names = cat[steps[0]].keys()
    for var in var_names:
        data = []
        step_data = []
        for step in steps:
            hfile[str(step)+"/"+var] = cat[step][var]
    hfile.close()

def cat_hdf5(file_loc,nums,output):
    print file_loc
    print output
    print nums
    outfile= h5py.File(output,'w')
    infiles = []
    for i in range(0,nums):
        infiles.append(h5py.File(file_loc.replace("${num}",str(i)),'r') )
    steps = infiles[0].keys()
    for step in steps:
        print "\n\n================"
        print "step: ",step
        keys = infiles[0][str(step)].keys()
        for key in keys:
            print "key: ",key
            data = []
            for i in range(0,nums):
                print "file num:",i
                data.append(infiles[i][step][key].value)
            data_tot = np.concatenate(data)
            var_path = "/"+str(step)+"/"+key
            print var_path
            outfile[var_path]=data_tot
            outfile.flush()
        for i in range(0,nums):
            infiles[i].flush()

pos_cat = dtk.Catalog(gio_loc)
pos_cat.add_steps(gltcs_steps)
var_names = ['x','y','z','vx','vy','vz',
             'nodeIndex','infall_halo_mass','infall_halo_tag',
             'hostIndex','host_halo_mass','host_halo_tag',
             'host_halo_x','host_halo_y','host_halo_z',
             'positionType']
var_names2 = ['x','y','z','vx','vy','vz',
              'nodeIndex','infallHaloMass','infallHaloTag',
              'hostIndex','hostHaloMass','hostHaloTag',
              'hostHaloX','hostHaloY','hostHaloZ',
              'positionType']
for i in range(0,len(var_names)):
    pos_cat.add_var(var_names[i],as_name=var_names2[i])

print "reading positions...",
start_read = time.time()    
#pos_cat.read_gio() TODO: remove this
print "done.", time.time()-start_read

gal_cat = dtk.Catalog(gltcs_loc)
gal_cat.add_steps(gltcs_steps,in_file_steps=gltcs_internal)
gal_cat.add_var_step_replace("${step_var_name}",gltcs_steps,gltcs_str_z)

#gal_cat.add_var("/Outputs/Output${step}/nodeData/nodeIndex",as_name='nodeIndex')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_u:rest:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_u:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_g:rest:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_g:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_r:rest:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_r:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_i:rest:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_i:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_z:rest:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_z:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_u:rest:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_u:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_g:rest:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_g:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_r:rest:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_r:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_i:rest:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_i:rest')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_z:rest:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_z:rest')

gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_u:observed:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_u:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_g:observed:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_g:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_r:observed:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_r:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_i:observed:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_i:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_z:observed:z${step_var_name}",as_name='spheroidLuminositiesStellar:SDSS_z:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_u:observed:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_u:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_g:observed:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_g:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_r:observed:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_r:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_i:observed:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_i:observed')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_z:observed:z${step_var_name}",as_name='diskLuminositiesStellar:SDSS_z:observed')

file_loc = gltcs_loc.replace("${subfile}",str(0))
hfile = h5py.File(file_loc,'r')
keys_tmp = hfile["/Outputs/Output1/nodeData/"].keys()
keys = []

#don't copy these table that have these words
avoids = ["LuminositiesStellar","hotHalo","basic","blackHoleCount","blackHoleJetPower","blackHoleRadiativeEfficiency","blackHoleSpin","darkMatterProfileScale","indicesBranchTip","position","satellite","siblingIndex","merge"]
for key in keys_tmp:
    if any(x in  key for x in avoids):
        print "throw",key
    else:
        print "keep ",key
        keys.append(key)
        gal_cat.add_var("/Outputs/Output${step}/nodeData/"+key,as_name=key)


keys.append('diskLuminositiesStellar:SDSS_u:rest')
keys.append('diskLuminositiesStellar:SDSS_g:rest')
keys.append('diskLuminositiesStellar:SDSS_r:rest')
keys.append('diskLuminositiesStellar:SDSS_i:rest')
keys.append('diskLuminositiesStellar:SDSS_z:rest')
keys.append('spheriodLuminositiesStellar:SDSS_u:rest')
keys.append('spheriodLuminositiesStellar:SDSS_g:rest')
keys.append('spheriodLuminositiesStellar:SDSS_r:rest')
keys.append('spheriodLuminositiesStellar:SDSS_i:rest')
keys.append('spheriodLuminositiesStellar:SDSS_z:rest')

keys.append('diskLuminositiesStellar:SDSS_u:observed')
keys.append('diskLuminositiesStellar:SDSS_g:observed')
keys.append('diskLuminositiesStellar:SDSS_r:observed')
keys.append('diskLuminositiesStellar:SDSS_i:observed')
keys.append('diskLuminositiesStellar:SDSS_z:observed')
keys.append('spheriodLuminositiesStellar:SDSS_u:observed')
keys.append('spheriodLuminositiesStellar:SDSS_g:observed')
keys.append('spheriodLuminositiesStellar:SDSS_r:observed')
keys.append('spheriodLuminositiesStellar:SDSS_i:observed')
keys.append('spheriodLuminositiesStellar:SDSS_z:observed')

final_create = True
results = []
tmp_file_num =0
del_list = ['magnitude:SDSS_u:rest','magnitude:SDSS_g:rest','magnitude:SDSS_r:rest','magnitude:SDSS_i:rest','magnitude:SDSS_z:rest', 
            'magnitude:SDSS_u:observed','magnitude:SDSS_g:observed','magnitude:SDSS_r:observed','magnitude:SDSS_i:observed','magnitude:SDSS_z:observed',
            'totalLuminositiesStellar:SDSS_u:observed','totalLuminositiesStellar:SDSS_g:observed','totalLuminositiesStellar:SDSS_r:observed','totalLuminositiesStellar:SDSS_i:observed','totalLuminositiesStellar:SDSS_z:observed',
            'totalLuminositiesStellar:SDSS_u:rest','totalLuminositiesStellar:SDSS_g:rest','totalLuminositiesStellar:SDSS_r:rest''totalLuminositiesStellar:SDSS_i:rest','totalLuminositiesStellar:SDSS_z:rest',
            'totalStarFormationRate','totalMassStellar','totalAbundancesStellarMetals','totalMassGas','totalAbundancesGasMetals','totalAngularMomentum',
            'totalLuminositiesStellar:SDSS_u:rest','totalLuminositiesStellar:SDSS_g:rest','totalLuminositiesStellar:SDSS_r:rest','totalLuminositiesStellar:SDSS_i:rest','totalLuminositiesStellar:SDSS_z:rest']

#for i in range(0,gltcs_file_num): TODO: remove
for i in range(0,0):
    start = time.time()
    gal_cat.add_subfile(i)
    if(i%gltcs_file_step!=gltcs_file_step-1 and not (i==gltcs_file_num-1)):
        #print "time: ",(time.time()-start)
        continue
    print "=============================================="
    print "working on file: ",i
    print "=============================================="
    start_read = time.time()
    gal_cat.delete_vars(del_list)
    gal_cat.read_hdf5_no_step_file()
    gal_cat.clear_subfiles()
    print "adding new vars"
    for step in gltcs_steps:
        print "step: ",step
        gal_cat[step]['magnitude:SDSS_u:rest']=-2.5*np.log10(gal_cat[step]['diskLuminositiesStellar:SDSS_u:rest']+gal_cat[step]['spheroidLuminositiesStellar:SDSS_u:rest'])
        gal_cat[step]['magnitude:SDSS_g:rest']=-2.5*np.log10(gal_cat[step]['diskLuminositiesStellar:SDSS_g:rest']+gal_cat[step]['spheroidLuminositiesStellar:SDSS_g:rest'])
        gal_cat[step]['magnitude:SDSS_r:rest']=-2.5*np.log10(gal_cat[step]['diskLuminositiesStellar:SDSS_r:rest']+gal_cat[step]['spheroidLuminositiesStellar:SDSS_r:rest'])
        gal_cat[step]['magnitude:SDSS_i:rest']=-2.5*np.log10(gal_cat[step]['diskLuminositiesStellar:SDSS_i:rest']+gal_cat[step]['spheroidLuminositiesStellar:SDSS_i:rest'])
        gal_cat[step]['magnitude:SDSS_z:rest']=-2.5*np.log10(gal_cat[step]['diskLuminositiesStellar:SDSS_z:rest']+gal_cat[step]['spheroidLuminositiesStellar:SDSS_z:rest'])
        for key in gal_cat[step].keys():
            if ("disk" in key) and ("Radius" not in key) and ("Velocity" not in key):
                sum_key = key.replace("disk","total")
                gal_cat[step][sum_key]=gal_cat[step][key]+gal_cat[step][key.replace("disk","spheroid")]

    gal_cat.refresh_vars()
    print "data shape: ",gal_cat[step][gal_cat[step].keys()[0]].size
   
            
    

    start_merge=time.time()
    print "mering one sided..."
    result = dtk.Catalog()
    result.quick_join(pos_cat,gal_cat,'nodeIndex')
    
    for step in gltcs_steps:
        redshift = stepz.get_z(step)
        dl = z_to_dl(redshift)
        adjust = -2.5*np.log10(1+redshift)+5*np.log10(dl)+25.0
        result[step]['magnitude:SDSS_u:observed']= adjust -2.5*np.log10(result[step]['diskLuminositiesStellar:SDSS_u:observed']+result[step]['spheroidLuminositiesStellar:SDSS_u:observed'])
        result[step]['magnitude:SDSS_g:observed']= adjust -2.5*np.log10(result[step]['diskLuminositiesStellar:SDSS_g:observed']+result[step]['spheroidLuminositiesStellar:SDSS_g:observed'])
        result[step]['magnitude:SDSS_r:observed']= adjust -2.5*np.log10(result[step]['diskLuminositiesStellar:SDSS_r:observed']+result[step]['spheroidLuminositiesStellar:SDSS_r:observed'])
        result[step]['magnitude:SDSS_i:observed']= adjust -2.5*np.log10(result[step]['diskLuminositiesStellar:SDSS_i:observed']+result[step]['spheroidLuminositiesStellar:SDSS_i:observed'])
        result[step]['magnitude:SDSS_z:observed']= adjust -2.5*np.log10(result[step]['diskLuminositiesStellar:SDSS_z:observed']+result[step]['spheroidLuminositiesStellar:SDSS_z:observed'])
    if(tmp_to_disk):
        print "writing to subfile to disk"
        catalog_to_hdf5(result,tmp_file_loc.replace("${num}",str(tmp_file_num)))
        tmp_file_num += 1
    else:
        results.append(result)

    end = time.time()
    print "time: ",(end-start),"\n\n\n"

print "merging everything together..."

if(tmp_to_disk):
    print "here we go to cat"
#    cat_hdf5(tmp_file_loc,tmp_file_num,output_box) TODO: remove
    cat_hdf5(tmp_file_loc,8,output_box)
else:
    total_result = pd.concat(results,ignore_index=True)
    print pos_df.shape, "vs", total_result.shape 
    print "diff", pos_df.shape[0]-total_result.shape[0]
    pandas_to_hdf5(total_result,output_box)

total_end = time.time()
print "done. Time:",(total_end - total_start)
