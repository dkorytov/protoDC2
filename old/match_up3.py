#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import sys
import time
import dtk
pd.set_option('display.expand_frame_repr', False)


param = dtk.Param(sys.argv[1])

lc_loc = param.get_string("lc_loc")
gltcs_loc = param.get_string("gltcs_loc")
gltcs_steps = param.get_int_list("gltcs_steps")
gltcs_internal = param.get_int_list("gltcs_internal")
gltcs_str_z    = param.get_string_list("gltcs_str_z")
gltcs_file_num = param.get_int("gltcs_file_num")
output         = param.get_string("output")
def load_bins(dic, loc,step,var_name,dtype,as_name=None):
    file_name = loc.replace("${step}",str(step)).replace("${var_name}",var_name)
    data = np.fromfile(file_name,dtype=dtype)
    if(as_name == None):
        dic[var_name]=data
    else:
        dic[as_name]=data
    print "from file: ",file_name
    print "\t loaded: ",data.size

def add_step(dic,step):
    size = dic[dic.keys()[0]].size
    data = np.ones(size,dtype="int")*step;
    dic["step"]=data
    print "addes steps: ",data.size

def pandas_to_hdf5(df,output):
    hfile = h5py.File(output,'w')
    columns = list(df.columns.values)
    for col in columns:
        hfile[col] = df[col]


dfs = []
for step in gltcs_steps:
    print "working on step: ", step
    lc_dic = {}
    load_bins(lc_dic,lc_loc,step,"id",dtype='i8',as_name='nodeIndex')
    load_bins(lc_dic,lc_loc,step,"phi",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"theta",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"redshift",dtype='f4')
    load_bins(lc_dic,lc_loc,step,'x',dtype='f4')
    load_bins(lc_dic,lc_loc,step,'y',dtype='f4')
    load_bins(lc_dic,lc_loc,step,'z',dtype='f4')
    load_bins(lc_dic,lc_loc,step,'vx',dtype='f4')
    load_bins(lc_dic,lc_loc,step,'vy',dtype='f4')
    load_bins(lc_dic,lc_loc,step,'vz',dtype='f4')
    add_step(lc_dic,step)
    df = pd.DataFrame(lc_dic)
    dfs.append(df)

lc_df = pd.concat(dfs,ignore_index=True)

print lc_df.head()
print lc_df.shape

print lc_df.describe()
exit()
gal_cat = dtk.Catalog(gltcs_loc)
gal_cat.add_steps(gltcs_steps,in_file_steps=gltcs_internal)
gal_cat.add_var_step_replace("${step_var_name}",gltcs_steps,gltcs_str_z)

gal_cat.add_var("/Outputs/Output${step}/nodeData/nodeIndex",as_name='nodeIndex')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_u:rest:z${step_var_name}",as_name='lum_sdss_u_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_g:rest:z${step_var_name}",as_name='lum_sdss_g_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_r:rest:z${step_var_name}",as_name='lum_sdss_r_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_i:rest:z${step_var_name}",as_name='lum_sdss_i_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidLuminositiesStellar:SDSS_z:rest:z${step_var_name}",as_name='lum_sdss_z_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidMassStellar",as_name='stellar_mass_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidMassGas",as_name='gas_mass_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/spheroidStarFormationRate",as_name='SFR_spheriod')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_u:rest:z${step_var_name}",as_name='lum_sdss_u_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_g:rest:z${step_var_name}",as_name='lum_sdss_g_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_r:rest:z${step_var_name}",as_name='lum_sdss_r_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_i:rest:z${step_var_name}",as_name='lum_sdss_i_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskLuminositiesStellar:SDSS_z:rest:z${step_var_name}",as_name='lum_sdss_z_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskMassStellar",as_name='stellar_mass_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskMassGas",as_name='gas_mass_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/diskStarFormationRate",as_name='SFR_disk')
gal_cat.add_var("/Outputs/Output${step}/nodeData/nodeIndex",as_name='hostIndex')
gal_cat.add_var("/Outputs/Output${step}/nodeData/nodeIsIsolated",as_name="isIsolated")


final_create = True
for i in range(0,gltcs_file_num):
    print "\n\n\n=============================================="
    print "working on file: ",i
    print "=============================================="
    start = time.time()
    gal_cat.clear_subfiles()
    gal_cat.add_subfile(i)
    gal_cat.read_hdf5()

    dfs = []
    for step in gltcs_steps:
        gal_dic = {}
        for key in gal_cat[step].keys():
            gal_dic[key] = gal_cat[step][key]
            gal_dic["step"]=np.ones(gal_dic[gal_dic.keys()[0]].size)*step
    dfs.append(pd.DataFrame(gal_dic))

    gal_df = pd.concat(dfs,ignore_index=True)
    gal_df['stellar_mass']= gal_df["stellar_mass_disk"] + gal_df["stellar_mass_spheriod"]
    gal_df['gas_mass']    = gal_df["gas_mass_disk"]     + gal_df["gas_mass_spheriod"]
    gal_df['SFR']         = gal_df['SFR_disk']          + gal_df['SFR_spheriod']
    gal_df['mag_sdss_u']=-2.5*np.log10(gal_df['lum_sdss_u_disk']+gal_df['lum_sdss_u_disk'])
    gal_df['mag_sdss_g']=-2.5*np.log10(gal_df['lum_sdss_g_disk']+gal_df['lum_sdss_g_disk'])
    gal_df['mag_sdss_r']=-2.5*np.log10(gal_df['lum_sdss_r_disk']+gal_df['lum_sdss_r_disk'])
    gal_df['mag_sdss_i']=-2.5*np.log10(gal_df['lum_sdss_i_disk']+gal_df['lum_sdss_i_disk'])
    gal_df['mag_sdss_z']=-2.5*np.log10(gal_df['lum_sdss_z_disk']+gal_df['lum_sdss_z_disk'])
    print "mering one sided..."
    if(final_create):
        final = lc_df.merge(gal_df,on=['step',"nodeIndex"],how="left").set_index(['step','nodeIndex'])
        final_create = False
    else:
        gal_df = gal_df.set_index(['step','nodeIndex'])
        print gal_df.head()
        print final.head()
        final.update(gal_df,raise_conflict=True)

    print final.head()
    print final.shape
    end = time.time()
    print "time: ",(end-start)

print "merging everything together..."
#total_df = pd.merge(gal_df,lc_df,on=['nodeIndex','step'])

print total_df.head()
print total_df.shape

pandas_to_hdf5(total_df,"test.hdf5")
