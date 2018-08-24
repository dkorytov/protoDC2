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
pd.set_option('display.expand_frame_repr', False)
h=0.702
total_start = time.time()
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
stepz = dtk.StepZ(200,0,500)
zs = np.linspace(0,1.5,1000)
z_to_dl = interp1d(zs,cosmo.luminosity_distance(zs))
def load_bins(dic, loc,step,var_name,dtype,as_name=None):
    file_name = loc.replace("${step}",str(step)).replace("${var_name}",var_name)
    data = np.fromfile(file_name,dtype=dtype)
    if(as_name == None):
        dic[var_name]=data
    else:
        dic[as_name]=data

def add_step(dic,step):
    size = dic[dic.keys()[0]].size
    data = np.ones(size,dtype="int")*step;
    dic["step"]=data
    print "\tadded steps. data size:",data.size

def pandas_to_hdf5(df,output):
    hfile = h5py.File(output,'w')
    columns = list(df.columns.values)
    print("pandas to hdf5")
    for col in columns:
        print("col",col)
        hfile[col] = df[col][:]
    hfile.flush()
    hfile.close()

def pandas_from_hdf5(output):
    hfile = h5py.File(output,'r')
    keys = hfile.keys()
    dic = {}
    for key in keys:
        dic[key]=hfile[key].value
    return pd.DataFrame(dic)

def load_shear(file_loc,step):
    dic = {}
    file_loc = file_loc.replace("${step}",str(step))
    #    data = np.loadtxt(file_loc,skiprows=1,dtype='i8,f4,f4,f4,f4,f4,f4,f4')
    nodeIndex,tx1,tx2,zs,s1,s2,k0,m0 = np.loadtxt(file_loc,skiprows=1,unpack=True,dtype={'names': ('nodeIndex', 'tx1', 'tx2','zs','s1','s2','k0','m0'), 'formats': ('i8',       'f4',  'f4','f4','f4','f4','f4','f4')})
    dic['nodeIndex']=nodeIndex
    dic['shear1']=s1
    dic['shear2']=s2
    dic['k0']    =k0
    dic['m0']    =m0
    dic['step']  =np.ones(nodeIndex.size,dtype='i4')*step
    print "loaded: ", nodeIndex.size
    return dic

def load_gio(file_loc,step):
    file_step_loc = file_loc.replace("${step}",str(step))
    dic = {}
    dic['nodeIndex']   = dtk.gio_read(file_step_loc,'nodeIndex')
    dic['hostIndex']   = dtk.gio_read(file_step_loc,'hostIndex')
    dic['hostHaloTag'] = dtk.gio_read(file_step_loc,"host_halo_tag")
    dic['hostHaloMass']= dtk.gio_read(file_step_loc,'host_halo_mass')/h
    dic['infallHaloTag']=dtk.gio_read(file_step_loc,'infall_halo_tag')
    dic['infallHaloMass']=dtk.gio_read(file_step_loc,'infall_halo_mass')/h
    dic['step']=np.ones(dic['nodeIndex'].size,dtype='i4')*step
    dic['placementType']=dtk.gio_read(file_step_loc,'positionType')
    return dic

def cat_hdf5(file_loc,nums,output):
    keys = h5py.File(file_loc.replace("${num}","0"),'r').keys()
    hfile= h5py.File(output,'w')
    for key in keys:
        data = []
        for i in range(0,nums):
            htmp = h5py.File(file_loc.replace("${num}",str(i)),'r')
            data.append(htmp[key][:])
        data_tot = np.concatenate(data)
        hfile[key]=data_tot
        hfile.flush()

dfs = []
for step in gltcs_steps:
    print "working on step: ", step
    lc_dic = {}
    load_bins(lc_dic,lc_loc,step,"id",dtype='i8',as_name='nodeIndex')
    load_bins(lc_dic,lc_loc,step,"phi",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"theta",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"redshift",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"x",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"y",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"z",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"vx",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"vy",dtype='f4')
    load_bins(lc_dic,lc_loc,step,"vz",dtype='f4')

    add_step(lc_dic,step)
    df1 = pd.DataFrame(lc_dic)
    gio_df = pd.DataFrame(load_gio(gio_loc,step))
    print "\t",df1.shape,gio_df.shape, "merged into",
    df = pd.merge(df1,gio_df,on=['nodeIndex','step'])
    print df.shape
    dfs.append(df)

lc_df = pd.concat(dfs,ignore_index=True)
print lc_df.head()
print lc_df.shape
if(use_shear):
    print "loading shears"
    if(True):
        sh_dfs = []
        for step in gltcs_steps:
            print "working on step: ", step
            sh_dic = load_shear(shear_loc,step)
            sh_dfs.append(pd.DataFrame(sh_dic))

        sh_df = pd.concat(sh_dfs,ignore_index=True)
        pandas_to_hdf5(sh_df,"tmp/shear.hdf5")
    else:
        sh_df = pandas_from_hdf5("tmp/shear.hdf5")
    print "pre merger shapes of lc and sh"
    print lc_df.shape
    print sh_df.shape
    lc_df = pd.merge(lc_df,sh_df,on=['step','nodeIndex'])
    print "post merger shape: "
    print lc_df.shape
    print "done loading"

gal_cat = dtk.Catalog(gltcs_loc)
gal_cat.add_steps(gltcs_steps,in_file_steps=gltcs_internal)
gal_cat.add_var_step_replace("${step_var_name}",gltcs_steps,gltcs_str_z)


file_loc = gltcs_loc.replace("${subfile}",str(0))
hfile = h5py.File(file_loc,'r')
keys_tmp = hfile["/Outputs/Output32/nodeData/"].keys()
keys = []

#don't copy these table that have these words
avoids = ["hotHalo","basic","blackHoleCount","blackHoleJetPower","blackHoleRadiativeEfficiency","blackHoleSpin","darkMatterProfileScale","indicesBranchTip","position","satellite","siblingIndex","merge"]

for key in keys_tmp:
    if any(x in  key for x in avoids):
        print "throw",key
    else:
        if "LuminositiesStellar" in key:
            key = key[:-8] #remove the ":z0.0000" and replace with a string that will be changed
            #to the correct value at each step
            gal_cat.add_var("/Outputs/Output${step}/nodeData/"+key+":z${step_var_name}",as_name=key)
        else:
            print "keep ",key
            gal_cat.add_var("/Outputs/Output${step}/nodeData/"+key,as_name=key)
        keys.append(key)


final_create = True
results = []
tmp_file_num =0
print gltcs_file_num
for i in range(0,gltcs_file_num):
    start = time.time()
    gal_cat.add_subfile(i)
    if(i%gltcs_file_step!=gltcs_file_step-1 and not (i==gltcs_file_num-1)):
        #print "time: ",(time.time()-start)
        continue
    print "\n\n\n=============================================="
    print "working on file: ",i
    print "=============================================="
    start_read = time.time()
    gal_cat.read_hdf5_no_step_file()
    gal_cat.clear_subfiles()
    dfs = []
    for step in gltcs_steps:
        gal_dic = {}
        for key in gal_cat[step].keys():
            gal_dic[key] = gal_cat[step][key]
        gal_dic["step"]=np.ones(gal_dic[gal_dic.keys()[0]].size)*step
        dfs.append(pd.DataFrame(gal_dic))

    gal_df = pd.concat(dfs,ignore_index=True)
    gal_df['magnitude:SDSS_u:rest']=-2.5*np.log10(gal_df['diskLuminositiesStellar:SDSS_u:rest']+gal_df['spheroidLuminositiesStellar:SDSS_u:rest'])
    gal_df['magnitude:SDSS_g:rest']=-2.5*np.log10(gal_df['diskLuminositiesStellar:SDSS_g:rest']+gal_df['spheroidLuminositiesStellar:SDSS_g:rest'])
    gal_df['magnitude:SDSS_r:rest']=-2.5*np.log10(gal_df['diskLuminositiesStellar:SDSS_r:rest']+gal_df['spheroidLuminositiesStellar:SDSS_r:rest'])
    gal_df['magnitude:SDSS_i:rest']=-2.5*np.log10(gal_df['diskLuminositiesStellar:SDSS_i:rest']+gal_df['spheroidLuminositiesStellar:SDSS_i:rest'])
    gal_df['magnitude:SDSS_z:rest']=-2.5*np.log10(gal_df['diskLuminositiesStellar:SDSS_z:rest']+gal_df['spheroidLuminositiesStellar:SDSS_z:rest'])
  
    print "read in time: ",(time.time()-start_read)
    if(use_mr_cut):
        start_cut = time.time()
        slct = gal_df['magnitude:SDSS_r:rest']<mr_cut
        print "cutting on mr: ",np.float(slct.sum())/np.float(slct.size),slct.sum(),slct.size
        gal_df = gal_df[slct]
        print "cut time: ",(time.time()-start_cut)
    print "data shape: ",gal_df.shape
    for key in keys:
        if ("disk" in key) and ("Radius" not in key) and ("Velocity" not in key):
            sum_key = key.replace("disk","total")
            gal_df[sum_key]=gal_df[key]+gal_df[key.replace("disk","spheroid")]
            
    

    start_merge=time.time()
    print "merging one sided..."
    result = pd.merge(lc_df,gal_df,on=["step","nodeIndex"])
    
    redshifts = result['redshift']
    print "min redshift: ",np.min(redshifts)
    print "max redshift: ",np.max(redshifts)
    dl = z_to_dl(redshifts)
    print "done getting the luminosity distance..."
    adjust = -2.5*np.log10(1+redshifts)+5*np.log10(dl)+25.0
    result['magnitude:SDSS_u:observed']= adjust -2.5*np.log10(result['diskLuminositiesStellar:SDSS_u:observed']+result['spheroidLuminositiesStellar:SDSS_u:observed'])
    result['magnitude:SDSS_g:observed']= adjust -2.5*np.log10(result['diskLuminositiesStellar:SDSS_g:observed']+result['spheroidLuminositiesStellar:SDSS_g:observed'])
    result['magnitude:SDSS_r:observed']= adjust -2.5*np.log10(result['diskLuminositiesStellar:SDSS_r:observed']+result['spheroidLuminositiesStellar:SDSS_r:observed'])
    result['magnitude:SDSS_i:observed']= adjust -2.5*np.log10(result['diskLuminositiesStellar:SDSS_i:observed']+result['spheroidLuminositiesStellar:SDSS_i:observed'])
    result['magnitude:SDSS_z:observed']= adjust -2.5*np.log10(result['diskLuminositiesStellar:SDSS_z:observed']+result['spheroidLuminositiesStellar:SDSS_z:observed'])
   
    if(tmp_to_disk):
        print "writing to subfile to disk"
        pandas_to_hdf5(result,tmp_file_loc.replace("${num}",str(tmp_file_num)))
        tmp_file_num += 1
    else:
        results.append(result)

    end = time.time()
    print "time: ",(end-start)

print "merging everything together..."

if(tmp_to_disk):
    print "here we go to cat"
    cat_hdf5(tmp_file_loc,tmp_file_num,output)
else:
    total_result = pd.concat(results,ignore_index=True)
    print lc_df.shape, "vs", total_result.shape 
    print "diff", lc_df.shape[0]-total_result.shape[0]
    pandas_to_hdf5(total_result,output)

total_end = time.time()
print "done. Time:",(total_end - total_start)

import subprocess 
subprocess.Popen("./add_metadata.py "+output, shell=True)

